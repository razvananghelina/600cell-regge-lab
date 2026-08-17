#!/usr/bin/env python3
"""Target-blind direct Regge acceleration on projected 600-cell refinements.

Prior-art commit: 81db1ec.
Protocol commit: 23157e2; unit-normalization correction: 05d5685;
quadratic-lapse correction after the preserved first failure: e0fcda4.

The refined coefficients are produced without comparing them with the
continuum value.  A separate post-commit verifier performs that comparison.
"""

from collections import Counter
from itertools import combinations
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "gravity_600cell_projected_refinement_acceleration_blind.json"
PRIOR_ART_COMMIT = "81db1ec"
PROTOCOL_COMMIT = "23157e2"
NORMALIZATION_CORRECTION_COMMIT = "05d5685"
LAPSE_CORRECTION_COMMIT = "e0fcda4"
VARIANTS = (
    "legacy_float_shortest",
    "first_tie_rank_0",
    "first_tie_rank_1",
    "first_tie_rank_2",
)
ETAS_G = (0.04, 0.02, 0.01, 0.005)
ETAS_F = ETAS_G[:3]
A_SENTINELS = (0.0, -1.0, -2.0, -3.0)
PRIMARY_DERIVATIVE_STEP = 2e-5
SECONDARY_DERIVATIVE_STEP = 1e-5
PRIMARY_LAPSE_DERIVATIVE_STEP = 2e-3
SECONDARY_LAPSE_DERIVATIVE_STEP = 1e-3
EXACT_COARSE_RADIUS_COEFFICIENT = -0.5394897340206755
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


def coarse_tetrahedra(adjacency):
    neighbors = [set(np.flatnonzero(adjacency[index] > 0.5))
                 for index in range(len(adjacency))]
    result = []
    for i in range(len(adjacency)):
        for j in sorted(vertex for vertex in neighbors[i] if vertex > i):
            common_ij = neighbors[i] & neighbors[j]
            for k in sorted(vertex for vertex in common_ij if vertex > j):
                common_ijk = common_ij & neighbors[k]
                for ell in sorted(
                        vertex for vertex in common_ijk if vertex > k):
                    result.append((i, j, k, ell))
    return np.asarray(result, dtype=np.int64)


def mesh_edges(tetrahedra):
    return sorted({tuple(sorted(edge)) for tet in tetrahedra
                   for edge in combinations(tet, 2)})


def mesh_faces_with_counts(tetrahedra):
    return Counter(tuple(sorted(face)) for tet in tetrahedra
                   for face in combinations(tet, 3))


def candidate_record(new_positions, a, b, c, d, midpoint_index):
    def midpoint(left, right):
        return midpoint_index[tuple(sorted((left, right)))]

    ab, ac, ad = midpoint(a, b), midpoint(a, c), midpoint(a, d)
    bc, bd, cd = midpoint(b, c), midpoint(b, d), midpoint(c, d)
    candidates = (
        (ab, cd, (ac, ad, bd, bc)),
        (ac, bd, (ab, ad, cd, bc)),
        (ad, bc, (ab, ac, cd, bd)),
    )
    records = []
    for left, right, cycle in candidates:
        pair = (min(left, right), max(left, right))
        length = float(np.linalg.norm(
            new_positions[left]-new_positions[right]))
        records.append({
            "length": length,
            "pair": pair,
            "left": left,
            "right": right,
            "cycle": cycle,
        })
    return records, (ab, ac, ad, bc, bd, cd)


def red_refine(positions, tetrahedra, first_tie_rank=None):
    new_positions = [point.copy() for point in positions]
    midpoint_index = {}
    for edge in mesh_edges(tetrahedra):
        midpoint = positions[edge[0]]+positions[edge[1]]
        midpoint /= np.linalg.norm(midpoint)
        midpoint_index[edge] = len(new_positions)
        new_positions.append(midpoint)
    new_positions = np.asarray(new_positions)

    refined = []
    selected_lengths = []
    candidate_spreads = []
    selected_pairs = []
    for tetrahedron in tetrahedra:
        a, b, c, d = sorted(int(value) for value in tetrahedron)
        records, mids = candidate_record(
            new_positions, a, b, c, d, midpoint_index)
        ab, ac, ad, bc, bd, cd = mids
        refined.extend((
            (a, ab, ac, ad),
            (b, ab, bc, bd),
            (c, ac, bc, cd),
            (d, ad, bd, cd),
        ))

        if first_tie_rank is None:
            selected = min(records, key=lambda record: (
                record["length"], record["pair"]))
        else:
            selected = sorted(records, key=lambda record: record["pair"])[
                first_tie_rank]
        selected_lengths.append(selected["length"])
        candidate_spreads.append(
            max(item["length"] for item in records)
            - min(item["length"] for item in records))
        selected_pairs.append(selected["pair"])
        left, right, cycle = (
            selected["left"], selected["right"], selected["cycle"])
        for index in range(4):
            refined.append(
                (left, right, cycle[index], cycle[(index+1) % 4]))

    return (
        new_positions,
        np.asarray(refined, dtype=np.int64),
        {
            "selected_length_min_mean_max": [
                float(np.min(selected_lengths)),
                float(np.mean(selected_lengths)),
                float(np.max(selected_lengths)),
            ],
            "candidate_spread_min_mean_max": [
                float(np.min(candidate_spreads)),
                float(np.mean(candidate_spreads)),
                float(np.max(candidate_spreads)),
            ],
            "selected_pair_digest": hashlib.sha256(
                np.asarray(selected_pairs, dtype=np.int64).tobytes()
            ).hexdigest(),
        },
    )


def triangle_area(points):
    u = points[..., 1, :]-points[..., 0, :]
    v = points[..., 2, :]-points[..., 0, :]
    square = (
        np.sum(u*u, axis=-1)*np.sum(v*v, axis=-1)
        - np.sum(u*v, axis=-1)**2
    )
    return 0.5*np.sqrt(np.maximum(square, 0.0))


def tetra_volume(points):
    edges = points[..., 1:, :]-points[..., :1, :]
    gram = np.einsum("...ai,...bi->...ab", edges, edges)
    return np.sqrt(np.maximum(np.linalg.det(gram), 0.0))/6.0


def richardson_four(values):
    values = np.asarray(values, dtype=float)
    first = (4*values[1:]-values[:-1])/3
    second = (16*first[1:]-first[:-1])/15
    return first, second


def richardson_three(values):
    first, second = richardson_four(values)
    if len(second) != 1:
        raise ValueError("three-value Richardson input required")
    return first, float(second[0])


def affine_root(residuals):
    alpha = float(residuals[0.0])
    beta = alpha-float(residuals[-1.0])
    if beta == 0:
        return np.nan
    return -alpha/beta


def affine_residual(residuals):
    values = (
        float(residuals[-2.0])-2*float(residuals[-1.0])
        + float(residuals[0.0]),
        float(residuals[-3.0])-2*float(residuals[-2.0])
        + float(residuals[-1.0]),
    )
    scale = max(abs(float(item)) for item in residuals.values())
    return max(abs(item) for item in values)/max(scale, 1e-30)


def quadratic_dynamic_root(residuals):
    f_minus_one = float(residuals[-1.0])
    f_minus_two = float(residuals[-2.0])
    kappa_2 = (f_minus_two-2*f_minus_one)/2
    kappa_1 = kappa_2-f_minus_one
    root = np.nan if kappa_2 == 0 else -kappa_1/kappa_2
    return float(root), float(kappa_1), float(kappa_2)


def quadratic_residual(residuals):
    _, kappa_1, kappa_2 = quadratic_dynamic_root(residuals)
    predicted = 9*kappa_2-3*kappa_1
    observed = float(residuals[-3.0])
    scale = max(abs(predicted), abs(observed), 1e-30)
    return abs(observed-predicted)/scale


def static_lapse_residual(residuals):
    scale = max(abs(float(item)) for key, item in residuals.items() if key != 0)
    return abs(float(residuals[0.0]))/max(scale, 1e-30)


class CellularMesh:
    """Precomputed geometry and direct homothetic cellular action."""

    local_pairs = tuple(combinations(range(4), 2))

    def __init__(self, label, positions, tetrahedra):
        self.label = label
        positions = np.asarray(positions, dtype=float)
        self.positions = positions/np.linalg.norm(positions, axis=1)[:, None]
        self.tetrahedra = np.asarray(tetrahedra, dtype=np.int64)
        self.points = self.positions[self.tetrahedra]
        self.edges = np.asarray(mesh_edges(self.tetrahedra), dtype=np.int64)
        face_counts = mesh_faces_with_counts(self.tetrahedra)
        self.face_counts = face_counts
        self.faces = np.asarray(sorted(face_counts), dtype=np.int64)

        self.edge_lengths = np.linalg.norm(
            self.positions[self.edges[:, 0]]-
            self.positions[self.edges[:, 1]], axis=1)
        self.face_areas = triangle_area(self.positions[self.faces])
        self.tetra_volumes = tetra_volume(self.points)
        self.volume_bar = float(np.sum(self.tetra_volumes))
        self.s0 = float((2*np.pi**2/self.volume_bar)**(1/3))

        self.differences = self.points[:, :3, :]-self.points[:, 3:4, :]
        self.spatial_gram = np.einsum(
            "nai,nbi->nab", self.differences, self.differences)
        self.cross_shape = np.einsum(
            "nai,ni->na", self.differences, self.points[:, 3, :])
        self.local_edge_lengths = np.stack([
            np.linalg.norm(self.points[:, i]-self.points[:, j], axis=1)
            for i, j in self.local_pairs
        ], axis=1)
        self.local_face_areas = np.stack([
            triangle_area(self.points[:, [j for j in range(4) if j != i], :])
            for i in range(4)
        ], axis=1)

        self.spatial_angles = self._spatial_dihedral_angles()
        curvature_bar = (
            2*np.pi*np.sum(self.edge_lengths)
            - np.sum(self.local_edge_lengths*self.spatial_angles)
        )
        self.curvature = float(self.s0*curvature_bar)
        self.mass = float(self.curvature/(8*np.pi))
        self.maximum_imaginary_action = 0.0
        self.minimum_lorentzian_positive = np.inf
        self.maximum_lorentzian_negative = -np.inf

    def topology(self):
        return {
            "vertices": int(len(self.positions)),
            "edges": int(len(self.edges)),
            "faces": int(len(self.faces)),
            "tetrahedra": int(len(self.tetrahedra)),
            "euler_characteristic": int(
                len(self.positions)-len(self.edges)+len(self.faces)
                - len(self.tetrahedra)),
            "face_incidence_values": sorted(set(self.face_counts.values())),
        }

    def digest(self):
        hasher = hashlib.sha256()
        hasher.update(np.round(self.positions, 14).tobytes())
        canonical = np.sort(self.tetrahedra, axis=1)
        canonical = canonical[np.lexsort(canonical.T[::-1])]
        hasher.update(canonical.tobytes())
        return hasher.hexdigest()

    def _spatial_dihedral_angles(self):
        inverse = np.linalg.inv(self.spatial_gram)
        normals = np.zeros((len(self.tetrahedra), 4, 3), dtype=float)
        normals[:, 0, 0] = 1
        normals[:, 1, 1] = 1
        normals[:, 2, 2] = 1
        normals[:, 3, :] = -1

        def inner(left, right):
            return np.einsum("nki,nij,nkj->nk", left, inverse, right)

        result = []
        for i, j in self.local_pairs:
            k, ell = [value for value in range(4) if value not in (i, j)]
            left = normals[:, k:k+1, :]
            right = normals[:, ell:ell+1, :]
            cosine = -inner(left, right)/np.sqrt(
                inner(left, left)*inner(right, right))
            result.append(np.arccos(np.clip(cosine[:, 0], -1, 1)))
        return np.stack(result, axis=1)

    def _metric_and_normals(self, s_minus, s_plus, rho):
        dtype = np.result_type(s_minus, s_plus, rho)
        metric = np.empty((len(self.tetrahedra), 4, 4), dtype=dtype)
        metric[:, :3, :3] = s_minus*s_minus*self.spatial_gram
        cross = s_minus*(s_plus-s_minus)*self.cross_shape
        metric[:, :3, 3] = cross
        metric[:, 3, :3] = cross
        metric[:, 3, 3] = -rho
        inverse = np.linalg.inv(metric)

        normals = np.zeros((len(self.tetrahedra), 4, 4), dtype=dtype)
        normals[:, 0, 0] = 1
        normals[:, 1, 1] = 1
        normals[:, 2, 2] = 1
        normals[:, 3, :3] = -1
        normals[:, 3, 3] = s_plus/s_minus-1
        bottom = np.zeros((len(self.tetrahedra), 1, 4), dtype=dtype)
        bottom[:, :, 3] = 1
        return metric, inverse, normals, bottom, -bottom

    @staticmethod
    def _inner(left, inverse, right):
        return np.einsum("nki,nij,nkj->nk", left, inverse, right)

    def _angle(self, left, inverse, right, boundary=False):
        cross = self._inner(left, inverse, right)
        product = (
            self._inner(left, inverse, left)
            * self._inner(right, inverse, right)
        )
        if boundary:
            denominator = -1j*np.sqrt(-product.astype(complex))
        else:
            denominator = np.sqrt(product.astype(complex))
        return np.arccos(-cross/denominator), product

    def gravitational_action(self, s_minus, s_plus, rho, audit=False):
        metric, inverse, normals, bottom, top = self._metric_and_normals(
            s_minus, s_plus, rho)
        lateral_angles = []
        lateral_products = []
        for i, j in self.local_pairs:
            k, ell = [value for value in range(4) if value not in (i, j)]
            angle, product = self._angle(
                normals[:, k:k+1, :], inverse,
                normals[:, ell:ell+1, :], False)
            lateral_angles.append(angle[:, 0])
            lateral_products.append(product[:, 0])
        lateral_angles = np.stack(lateral_angles, axis=1)
        lateral_products = np.stack(lateral_products, axis=1)

        bottom_angles = []
        top_angles = []
        boundary_products = []
        for index in range(4):
            angle, product = self._angle(
                bottom, inverse, normals[:, index:index+1, :], True)
            bottom_angles.append(angle[:, 0])
            boundary_products.append(product[:, 0])
            angle, product = self._angle(
                top, inverse, normals[:, index:index+1, :], True)
            top_angles.append(angle[:, 0])
            boundary_products.append(product[:, 0])
        bottom_angles = np.stack(bottom_angles, axis=1)
        top_angles = np.stack(top_angles, axis=1)
        boundary_products = np.stack(boundary_products, axis=1)

        delta = s_plus-s_minus
        edge_areas = (
            1j*(s_minus+s_plus)/2*self.edge_lengths
            * np.sqrt(rho+delta*delta*self.edge_lengths**2/4)
        )
        local_edge_areas = (
            1j*(s_minus+s_plus)/2*self.local_edge_lengths
            * np.sqrt(rho+delta*delta*self.local_edge_lengths**2/4)
        )
        sum_hinges = (
            2*np.pi*np.sum(edge_areas)
            - np.sum(local_edge_areas*lateral_angles)
            + np.pi*s_minus*s_minus*np.sum(self.face_areas)
            - np.sum(s_minus*s_minus*self.local_face_areas*bottom_angles)
            + np.pi*s_plus*s_plus*np.sum(self.face_areas)
            - np.sum(s_plus*s_plus*self.local_face_areas*top_angles)
        )
        result = -1j*sum_hinges

        diagnostics = None
        if audit:
            if np.iscomplexobj(metric):
                raise ValueError("real metric required for inertia audit")
            eigenvalues = np.linalg.eigvalsh(metric)
            positive = eigenvalues[:, 1]
            negative = eigenvalues[:, 0]
            self.minimum_lorentzian_positive = min(
                self.minimum_lorentzian_positive, float(np.min(positive)))
            self.maximum_lorentzian_negative = max(
                self.maximum_lorentzian_negative, float(np.max(negative)))
            diagnostics = {
                "maximum_bottom_pi_over_2_residual": float(np.max(
                    np.abs(bottom_angles-np.pi/2))),
                "maximum_top_pi_over_2_residual": float(np.max(
                    np.abs(top_angles-np.pi/2))),
                "minimum_lateral_norm_product": float(np.min(
                    np.real(lateral_products))),
                "maximum_boundary_norm_product": float(np.max(
                    np.real(boundary_products))),
                "minimum_second_eigenvalue": float(np.min(positive)),
                "maximum_first_eigenvalue": float(np.max(negative)),
                "lorentzian_inertia_3_1": bool(
                    np.all(eigenvalues[:, 0] < 0)
                    and np.all(eigenvalues[:, 1:] > 0)),
            }
        if not any(np.iscomplexobj(value) for value in (
                s_minus, s_plus, rho)):
            self.maximum_imaginary_action = max(
                self.maximum_imaginary_action, float(abs(np.imag(result))))
        return result, diagnostics

    def total_action(self, s_minus, s_plus, rho):
        gravitational, _ = self.gravitational_action(
            s_minus, s_plus, rho, audit=False)
        return gravitational-8*np.pi*self.mass*np.sqrt(rho)

    def log_derivative(self, s_minus, s_plus, rho, coordinate, base_step):
        values = [s_minus, s_plus, rho]

        def centered(step):
            plus = [complex(value) for value in values]
            minus = [complex(value) for value in values]
            plus[coordinate] *= np.exp(1j*step)
            minus[coordinate] *= np.exp(-1j*step)
            return (
                self.total_action(*plus)-self.total_action(*minus)
            )/(2j*step)

        coarse = centered(base_step)
        fine = centered(base_step/2)
        return (4*fine-coarse)/3

    def real_log_derivative(self, s_minus, s_plus, rho, coordinate,
                            base_step):
        values = [float(s_minus), float(s_plus), float(rho)]

        def centered(step):
            plus = values.copy()
            minus = values.copy()
            plus[coordinate] *= np.exp(step)
            minus[coordinate] *= np.exp(-step)
            return (
                self.total_action(*plus)-self.total_action(*minus)
            )/(2*step)

        coarse = centered(base_step)
        fine = centered(base_step/2)
        return (4*fine-coarse)/3

    def seam_residual_table(self, derivative_step):
        seam_table = {a: [] for a in A_SENTINELS}
        maximum_derivative_imaginary = 0.0
        real_state_ok = True
        for eta in ETAS_G:
            rho = eta*eta
            previous_plus = 0.5*self.log_derivative(
                self.s0, self.s0, rho, 1, derivative_step)
            for a in A_SENTINELS:
                s_plus = self.s0*np.exp(a*eta*eta)
                _, state_audit = self.gravitational_action(
                    self.s0, s_plus, rho, audit=True)
                real_state_ok &= state_audit["lorentzian_inertia_3_1"]
                current_minus = 0.5*self.log_derivative(
                    self.s0, s_plus, rho, 0, derivative_step)
                seam = (previous_plus+current_minus)/eta
                seam_table[a].append(float(np.real(seam)))
                maximum_derivative_imaginary = max(
                    maximum_derivative_imaginary, float(abs(np.imag(seam))))
        return seam_table, {
            "maximum_scaled_derivative_imaginary": maximum_derivative_imaginary,
            "all_real_states_lorentzian": bool(real_state_ok),
        }

    def lapse_residual_table(self, derivative_step):
        lapse_table = {a: [] for a in A_SENTINELS}
        maximum_derivative_imaginary = 0.0
        for eta in ETAS_F:
            rho = eta*eta
            for a in A_SENTINELS:
                s_plus = self.s0*np.exp(a*eta*eta)
                lapse = self.real_log_derivative(
                    self.s0, s_plus, rho, 2, derivative_step)/eta**3
                lapse_table[a].append(float(np.real(lapse)))
                maximum_derivative_imaginary = max(
                    maximum_derivative_imaginary, float(abs(np.imag(lapse))))
        return lapse_table, {
            "maximum_scaled_derivative_imaginary": maximum_derivative_imaginary,
        }

    def coefficient_audit(self):
        seam, primary_diag = self.seam_residual_table(
            PRIMARY_DERIVATIVE_STEP)
        lapse, lapse_primary_diag = self.lapse_residual_table(
            PRIMARY_LAPSE_DERIVATIVE_STEP)
        seam_coarse = {}
        seam_fine = {}
        seam_truncations = {}
        for a in A_SENTINELS:
            _, second = richardson_four(seam[a])
            seam_coarse[a] = float(second[0])
            seam_fine[a] = float(second[1])
            seam_truncations[a] = float(abs(second[1]-second[0]))
        a_seam_coarse = affine_root(seam_coarse)
        a_seam = affine_root(seam_fine)

        lapse_limit = {}
        for a in A_SENTINELS:
            _, lapse_limit[a] = richardson_three(lapse[a])
        a_lapse, lapse_kappa_1, lapse_kappa_2 = quadratic_dynamic_root(
            lapse_limit)

        seam_secondary, secondary_diag = self.seam_residual_table(
            SECONDARY_DERIVATIVE_STEP)
        secondary_limit = {}
        for a in A_SENTINELS:
            _, second = richardson_four(seam_secondary[a])
            secondary_limit[a] = float(second[1])
        a_secondary = affine_root(secondary_limit)

        lapse_secondary, lapse_secondary_diag = self.lapse_residual_table(
            SECONDARY_LAPSE_DERIVATIVE_STEP)
        lapse_secondary_limit = {}
        for a in A_SENTINELS:
            _, lapse_secondary_limit[a] = richardson_three(
                lapse_secondary[a])
        a_lapse_secondary, _, _ = quadratic_dynamic_root(
            lapse_secondary_limit)

        return {
            "coefficient": float(a_seam),
            "coarse_extrapolation_coefficient": float(a_seam_coarse),
            "secondary_derivative_step_coefficient": float(a_secondary),
            "lapse_coefficient": float(a_lapse),
            "secondary_lapse_derivative_step_coefficient": float(
                a_lapse_secondary),
            "coefficient_truncation_difference": float(
                abs(a_seam-a_seam_coarse)),
            "derivative_step_difference": float(abs(a_seam-a_secondary)),
            "lapse_seam_difference": float(abs(a_seam-a_lapse)),
            "lapse_derivative_step_difference": float(
                abs(a_lapse-a_lapse_secondary)),
            "seam_affine_relative_residual": float(
                affine_residual(seam_fine)),
            "lapse_quadratic_relative_residual": float(
                quadratic_residual(lapse_limit)),
            "lapse_static_relative_residual": float(
                static_lapse_residual(lapse_limit)),
            "lapse_kappa_1": float(lapse_kappa_1),
            "lapse_kappa_2": float(lapse_kappa_2),
            "seam_residuals": {
                str(a): [float(item) for item in seam[a]]
                for a in A_SENTINELS
            },
            "lapse_residuals": {
                str(a): [float(item) for item in lapse[a]]
                for a in A_SENTINELS
            },
            "seam_limit_coarse": {
                str(a): float(seam_coarse[a]) for a in A_SENTINELS
            },
            "seam_limit_fine": {
                str(a): float(seam_fine[a]) for a in A_SENTINELS
            },
            "lapse_limit": {
                str(a): float(lapse_limit[a]) for a in A_SENTINELS
            },
            "secondary_lapse_limit": {
                str(a): float(lapse_secondary_limit[a])
                for a in A_SENTINELS
            },
            "maximum_seam_residual_truncation": float(
                max(seam_truncations.values())),
            "primary_diagnostics": primary_diag,
            "secondary_diagnostics": secondary_diag,
            "lapse_primary_diagnostics": lapse_primary_diag,
            "lapse_secondary_diagnostics": lapse_secondary_diag,
        }


def closed_regular_action(l_minus, l_plus, rho):
    delta = l_plus-l_minus
    h = np.sqrt(rho+delta**2/4)
    cosine = (delta**2+2*rho)/(2*(delta**2+3*rho))
    boost = delta/np.sqrt(8*(delta**2+3*rho))
    epsilon = 2*np.pi-5*np.arccos(cosine)
    return (
        360*(l_minus+l_plus)*h*epsilon
        + 600*np.sqrt(3)*(l_minus**2-l_plus**2)*np.arcsinh(boost)
    )


def mesh_record(mesh, coefficient, static_audit, relabel_error,
                static_action_error, static_total_error,
                static_lapse_relative):
    return {
        "topology": mesh.topology(),
        "mesh_digest": mesh.digest(),
        "volume_bar": mesh.volume_bar,
        "scale_for_unit_volume_radius": mesh.s0,
        "spatial_regge_curvature": mesh.curvature,
        "selected_total_dust_mass": mesh.mass,
        "minimum_tetra_volume": float(np.min(mesh.tetra_volumes)),
        "static_audit": static_audit,
        "static_action_relative_error": static_action_error,
        "static_total_action_relative_error": static_total_error,
        "static_lapse_relative_residual": static_lapse_relative,
        "relabel_action_relative_error": relabel_error,
        "maximum_imaginary_action": mesh.maximum_imaginary_action,
        "minimum_lorentzian_positive_eigenvalue": (
            mesh.minimum_lorentzian_positive),
        "maximum_lorentzian_negative_eigenvalue": (
            mesh.maximum_lorentzian_negative),
        "coefficient_audit": coefficient,
    }


def evaluate_mesh(label, positions, tetrahedra):
    print(f"[INFO] evaluating {label}", flush=True)
    mesh = CellularMesh(label, positions, tetrahedra)
    eta_static = 0.02
    static_grav, static_audit = mesh.gravitational_action(
        mesh.s0, mesh.s0, eta_static**2, audit=True)
    expected_static = eta_static*mesh.curvature
    static_action_error = float(
        abs(static_grav-expected_static)/max(abs(expected_static), 1e-30))
    static_total = mesh.total_action(mesh.s0, mesh.s0, eta_static**2)
    static_total_error = float(
        abs(static_total)/max(abs(expected_static), 1e-30))
    static_lapse = mesh.log_derivative(
        mesh.s0, mesh.s0, eta_static**2, 2, PRIMARY_DERIVATIVE_STEP)
    static_lapse_relative = float(
        abs(static_lapse)/max(abs(expected_static), 1e-30))

    permutation = np.array((2, 0, 3, 1), dtype=np.int64)
    reordered = CellularMesh(
        label+"_relabel", mesh.positions, mesh.tetrahedra[:, permutation])
    nonstatic = (mesh.s0, mesh.s0*0.97, 0.02)
    original_action, _ = mesh.gravitational_action(*nonstatic, audit=True)
    reordered_action, _ = reordered.gravitational_action(*nonstatic, audit=True)
    relabel_error = float(
        abs(original_action-reordered_action)
        / max(abs(original_action), abs(reordered_action), 1e-30))
    del reordered

    coefficient = mesh.coefficient_audit()
    record = mesh_record(
        mesh, coefficient, static_audit, relabel_error,
        static_action_error, static_total_error, static_lapse_relative)
    internal_ok = bool(
        mesh.topology()["euler_characteristic"] == 0
        and mesh.topology()["face_incidence_values"] == [2]
        and np.min(mesh.tetra_volumes) > 0
        and static_audit["lorentzian_inertia_3_1"]
        and static_audit["maximum_bottom_pi_over_2_residual"] < 5e-9
        and static_audit["maximum_top_pi_over_2_residual"] < 5e-9
        and static_action_error < 5e-9
        and static_total_error < 5e-8
        and static_lapse_relative < 5e-8
        and relabel_error < 5e-9
        and coefficient["coefficient_truncation_difference"] < 2e-6
        and coefficient["derivative_step_difference"] < 2e-6
        and coefficient["lapse_derivative_step_difference"] < 5e-5
        and coefficient["lapse_seam_difference"] < 5e-5
        and coefficient["seam_affine_relative_residual"] < 2e-6
        and coefficient["lapse_quadratic_relative_residual"] < 2e-6
        and coefficient["lapse_static_relative_residual"] < 2e-6
        and coefficient["primary_diagnostics"]["all_real_states_lorentzian"]
        and coefficient["secondary_diagnostics"]["all_real_states_lorentzian"]
    )
    return record, internal_ok


print("="*78)
print("BLIND DIRECT REGGE ACCELERATION ON PROJECTED 600-CELL REFINEMENTS")
print("="*78)

protocol_ok = bool(
    PRIOR_ART_COMMIT == "81db1ec"
    and PROTOCOL_COMMIT == "23157e2"
    and NORMALIZATION_CORRECTION_COMMIT == "05d5685"
    and LAPSE_CORRECTION_COMMIT == "e0fcda4"
    and len(VARIANTS) == 4
)
check("the frozen prior-art and protocol commits are named", protocol_ok)

source_positions, adjacency, _ = build_600cell()
source_norm_residual = float(np.max(
    np.abs(np.linalg.norm(source_positions, axis=1)-1)))
positions0 = source_positions/np.linalg.norm(source_positions, axis=1)[:, None]
tetrahedra0 = coarse_tetrahedra(adjacency)

towers = {}
refinement_audits = {}
for variant in VARIANTS:
    rank = None if variant == "legacy_float_shortest" else int(variant[-1])
    positions1, tetrahedra1, audit1 = red_refine(
        positions0, tetrahedra0, first_tie_rank=rank)
    positions2, tetrahedra2, audit2 = red_refine(
        positions1, tetrahedra1, first_tie_rank=None)
    towers[variant] = (
        (positions0, tetrahedra0),
        (positions1, tetrahedra1),
        (positions2, tetrahedra2),
    )
    refinement_audits[variant] = {"level1": audit1, "level2": audit2}

expected_sizes = ((120, 600), (840, 4800), (6480, 38400))
size_ok = all(
    tuple((len(pos), len(tet)) for pos, tet in towers[variant])
    == expected_sizes for variant in VARIANTS)
check("all four registered towers have the frozen carrier sizes", size_ok)

level1_exact_tie_scale = max(
    refinement_audits[variant]["level1"]
    ["candidate_spread_min_mean_max"][2]
    for variant in VARIANTS)
pair_digests = {
    refinement_audits[variant]["level1"]["selected_pair_digest"]
    for variant in VARIANTS
}
diagonal_audit_ok = bool(
    level1_exact_tie_scale < 1e-8 and len(pair_digests) >= 3)
check(
    "the exact first-level diagonal tie is exposed rather than hidden",
    diagonal_audit_ok,
    f"max float spread={level1_exact_tie_scale:.3e}; "
    f"distinct pair digests={len(pair_digests)}",
)

# Evaluate the common coarse carrier once.
records = {variant: {} for variant in VARIANTS}
coarse_record, coarse_internal_ok = evaluate_mesh(
    "level0", positions0, tetrahedra0)
for variant in VARIANTS:
    records[variant]["level0"] = coarse_record

coarse_mesh = CellularMesh("level0_closed_control", positions0, tetrahedra0)
phi = (1+np.sqrt(5))/2
closed_errors = []
closed_imaginaries = []
for l_minus, l_plus, rho in (
        (1.0, 1.0, 0.25), (1.0, 0.8, 0.1), (1.0, 1.2, 0.1)):
    direct, _ = coarse_mesh.gravitational_action(
        phi*l_minus, phi*l_plus, rho, audit=True)
    expected = closed_regular_action(l_minus, l_plus, rho)
    closed_errors.append(float(
        abs(direct-expected)/max(abs(expected), 1e-30)))
    closed_imaginaries.append(float(abs(np.imag(direct))))
closed_action_ok = bool(
    max(closed_errors) < 5e-10 and max(closed_imaginaries) < 5e-9)
check(
    "the irregular action reproduces the exact regular closed action",
    closed_action_ok,
    f"max relative={max(closed_errors):.3e}; "
    f"max imaginary={max(closed_imaginaries):.3e}",
)

zeta = (np.pi**2*np.sqrt(2)/50)**(1/3)
epsilon3 = 2*np.pi-5*np.arccos(1/3)
scaled_edges = coarse_mesh.s0*coarse_mesh.edge_lengths
volume_radius_ok = bool(
    np.max(abs(scaled_edges-zeta)) < 5e-10
    and abs(coarse_mesh.curvature-720*zeta*epsilon3)
        / abs(720*zeta*epsilon3) < 5e-9)
check(
    "the volume radius and spatial curvature recover the exact coarse map",
    volume_radius_ok,
    f"edge spread={np.max(abs(scaled_edges-zeta)):.3e}; "
    f"C={coarse_mesh.curvature:.15g}",
)

coarse_coefficient_ok = bool(
    abs(coarse_record["coefficient_audit"]["coefficient"]
        - EXACT_COARSE_RADIUS_COEFFICIENT) < 5e-6)
check(
    "the frozen estimator recovers the exact coarse acceleration",
    coarse_coefficient_ok,
    f"a0={coarse_record['coefficient_audit']['coefficient']:.15g}",
)

all_internal_ok = coarse_internal_ok
for variant in VARIANTS:
    for level in (1, 2):
        positions, tetrahedra = towers[variant][level]
        record, ok = evaluate_mesh(
            f"{variant}/level{level}", positions, tetrahedra)
        records[variant][f"level{level}"] = record
        all_internal_ok &= ok
        print(
            f"[BLIND] {variant} level{level}: "
            f"a={record['coefficient_audit']['coefficient']:.15g}",
            flush=True,
        )

check(
    "every refined carrier passes topology, action, branch and estimator gates",
    all_internal_ok,
)

all_coefficients_finite = all(
    np.isfinite(records[variant][f"level{level}"]
                ["coefficient_audit"]["coefficient"])
    for variant in VARIANTS for level in range(3)
)
check("all blind coefficients are finite", all_coefficients_finite)

outcome = (
    "PROJECTED_REGGE_ACCELERATION_COEFFICIENTS_DERIVED"
    if all((
        protocol_ok, size_ok, diagonal_audit_ok, closed_action_ok,
        volume_radius_ok, coarse_coefficient_ok, all_internal_ok,
        all_coefficients_finite,
    ))
    else "PROJECTED_REGGE_ACCELERATION_COEFFICIENTS_OPEN"
)

payload = {
    "protocol": (
        "direct coefficient production; no refined continuum-target "
        "comparison performed"
    ),
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "normalization_correction_commit": NORMALIZATION_CORRECTION_COMMIT,
    "lapse_correction_commit": LAPSE_CORRECTION_COMMIT,
    "source_maximum_unit_norm_residual": source_norm_residual,
    "variants": list(VARIANTS),
    "eta_values_seam": list(ETAS_G),
    "eta_values_lapse": list(ETAS_F),
    "a_sentinels": list(A_SENTINELS),
    "primary_derivative_step": PRIMARY_DERIVATIVE_STEP,
    "secondary_derivative_step": SECONDARY_DERIVATIVE_STEP,
    "primary_lapse_derivative_step": PRIMARY_LAPSE_DERIVATIVE_STEP,
    "secondary_lapse_derivative_step": SECONDARY_LAPSE_DERIVATIVE_STEP,
    "refinement_audits": refinement_audits,
    "records": records,
    "closed_regular_control": {
        "maximum_relative_error": max(closed_errors),
        "maximum_imaginary_contamination": max(closed_imaginaries),
    },
    "checks": {
        "protocol_ok": protocol_ok,
        "size_ok": size_ok,
        "diagonal_audit_ok": diagonal_audit_ok,
        "closed_action_ok": closed_action_ok,
        "volume_radius_ok": volume_radius_ok,
        "coarse_coefficient_ok": coarse_coefficient_ok,
        "all_internal_ok": all_internal_ok,
        "all_coefficients_finite": all_coefficients_finite,
    },
    "tests": tests+1,
    "passed": passed+int(outcome.endswith("DERIVED")),
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

check(
    "the Stage-A outcome follows the frozen internal hierarchy",
    outcome == "PROJECTED_REGGE_ACCELERATION_COEFFICIENTS_DERIVED",
    f"outcome={outcome}",
)

print("\nBlind coefficients (no refined continuum comparison):")
print(f"  common level0: {coarse_record['coefficient_audit']['coefficient']:.15g}")
for variant in VARIANTS:
    print(
        f"  {variant}: "
        + ", ".join(
            f"level{level}="
            f"{records[variant][f'level{level}']['coefficient_audit']['coefficient']:.15g}"
            for level in (1, 2)
        )
    )
print(f"\nSummary: {passed}/{tests} checks passed")
print(f"Outcome: {outcome}")
print(f"Artifact: {OUTPUT}")
raise SystemExit(0 if passed == tests else 1)

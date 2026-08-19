#!/usr/bin/env python3
"""Orbit-compressed adversarial audit of canonical-carrier acceleration."""

import ast
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
import hashlib
import json
from pathlib import Path
import sys

import networkx as nx
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_600cell_projected_rank_edgewise_acceleration_adversarial_corrected.json"
PROTOCOL_COMMIT = "3b7bd6c"
PROTOCOL_CORRECTION_COMMIT = "1c3318d"
PRIMARY_BLIND = HERE / "gravity_600cell_projected_rank_edgewise_acceleration_blind.json"
PRIMARY_COMPARISON = HERE / "gravity_600cell_projected_rank_edgewise_acceleration_comparison.json"
CARRIER_ARTIFACT = HERE / "gravity_600cell_projected_rank_edgewise_carrier.json"
ACTION_SOURCE = HERE / "verify_gravity_600cell_projected_refinement_acceleration_blind.py"
INPUT_HASHES = {
    PRIMARY_BLIND:
        "2059620f22cfbd8eac8abe6f2c7536924128d37f47a430bf773e34a9aead93a2",
    PRIMARY_COMPARISON:
        "132a81fe03ee67dbe95b91a68910f9212db88f8c6104b23f7d3f3f422939f5a4",
    CARRIER_ARTIFACT:
        "b57955b85a972df00b5673ddf7ee295757848f5afb43314857cf3de2dc85ac84",
}
ETAS_SEAM = (0.05, 0.025, 0.0125, 0.00625)
ETAS_LAPSE = ETAS_SEAM[:3]
SENTINELS = (0.0, -0.75, -1.5, -2.25)
SEAM_STEPS = (4e-4, 2e-4)
LAPSE_STEPS = (4e-3, 2e-3)
REGULAR_CONTROL = -0.5394897340206755
LOCAL_PAIRS = tuple(combinations(range(4), 2))
ROUND_DIGITS = 13
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"[{'PASS' if condition else 'FAIL'}] {label}", flush=True)
    if detail:
        print(f"       {detail}", flush=True)
    return condition


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def all_simplices(top):
    return tuple(
        tuple(sorted({tuple(sorted(face)) for tetrahedron in top
                      for face in combinations(tetrahedron, degree+1)}))
        for degree in range(4)
    )


def networkx_tetrahedra(adjacency):
    graph = nx.from_numpy_array(adjacency)
    cliques = [tuple(sorted(clique)) for clique in nx.find_cliques(graph)]
    return tuple(sorted(cliques))


def projected_barycentric(vertices, coarse_cells, coarse_top):
    vertex_cells = tuple(cell for layer in coarse_cells for cell in layer)
    index = {cell: position for position, cell in enumerate(vertex_cells)}
    raw = np.asarray([vertices[list(cell)].mean(axis=0)
                      for cell in vertex_cells])
    positions = raw/np.linalg.norm(raw, axis=1)[:, None]
    top = []
    for tetrahedron in coarse_top:
        for ordering in permutations(tetrahedron):
            flags = (
                (ordering[0],), tuple(sorted(ordering[:2])),
                tuple(sorted(ordering[:3])), tetrahedron,
            )
            top.append(tuple(index[flag] for flag in flags))
    return positions, tuple(top)


DIRECT_CHILDREN = (
    ("v0", "m01", "m02", "m03"),
    ("v1", "m01", "m12", "m13"),
    ("v2", "m02", "m12", "m23"),
    ("v3", "m03", "m13", "m23"),
    ("m01", "m02", "m03", "m13"),
    ("m01", "m02", "m12", "m13"),
    ("m02", "m03", "m13", "m23"),
    ("m02", "m12", "m13", "m23"),
)


def direct_rank_split(base_positions, base_top):
    top_keys = []
    vertex_keys = set()
    for chamber in base_top:
        names = {f"v{rank}": (chamber[rank], chamber[rank])
                 for rank in range(4)}
        for left, right in LOCAL_PAIRS:
            names[f"m{left}{right}"] = tuple(sorted(
                (chamber[left], chamber[right])
            ))
        for child in DIRECT_CHILDREN:
            keys = tuple(sorted(names[name] for name in child))
            top_keys.append(keys)
            vertex_keys.update(keys)
    ordered = tuple(sorted(vertex_keys))
    index = {key: position for position, key in enumerate(ordered)}
    raw = np.asarray([base_positions[list(key)].mean(axis=0)
                      for key in ordered])
    positions = raw/np.linalg.norm(raw, axis=1)[:, None]
    top = tuple(tuple(sorted(index[key] for key in tetrahedron))
                for tetrahedron in top_keys)
    return positions, top


def squared_distances(points):
    differences = points[:, None, :]-points[None, :, :]
    return np.einsum("ija,ija->ij", differences, differences)


def edge_signature(points):
    raw = (float(np.sum((points[0]-points[1])**2)),)
    return tuple(round(value, ROUND_DIGITS) for value in raw), raw


def face_signature(points):
    matrix = squared_distances(points)
    raw = tuple(sorted(float(matrix[left, right])
                       for left, right in combinations(range(3), 2)))
    return tuple(round(value, ROUND_DIGITS) for value in raw), raw


def tetra_signature(points):
    matrix = squared_distances(points)
    candidates = []
    for permutation in permutations(range(4)):
        candidates.append(tuple(float(matrix[permutation[left],
                                             permutation[right]])
                                for left, right in LOCAL_PAIRS))
    raw = min(candidates)
    return tuple(round(value, ROUND_DIGITS) for value in raw), raw


def classify_entities(positions, entities, signature):
    classes = {}
    maximum_residual = 0.0
    for entity in entities:
        points = positions[list(entity)]
        key, raw = signature(points)
        if key not in classes:
            classes[key] = {
                "representative_entity": tuple(entity),
                "representative_points": points.copy(),
                "representative_raw": raw,
                "multiplicity": 0,
            }
        record = classes[key]
        record["multiplicity"] += 1
        maximum_residual = max(maximum_residual, max(
            abs(left-right) for left, right in zip(
                raw, record["representative_raw"]
            )
        ))
    return tuple(classes[key] for key in sorted(classes)), maximum_residual


def triangle_area(points):
    u = points[..., 1, :]-points[..., 0, :]
    v = points[..., 2, :]-points[..., 0, :]
    square = (np.sum(u*u, axis=-1)*np.sum(v*v, axis=-1)
              - np.sum(u*v, axis=-1)**2)
    return 0.5*np.sqrt(np.maximum(square, 0.0))


def tetra_volume(points):
    edges = points[..., 1:, :]-points[..., :1, :]
    gram = np.einsum("...ai,...bi->...ab", edges, edges)
    return np.sqrt(np.maximum(np.linalg.det(gram), 0.0))/6.0


class ReducedMesh:
    def __init__(self, label, positions, top):
        self.label = label
        self.positions = np.asarray(positions, dtype=float)
        self.top = tuple(tuple(tetrahedron) for tetrahedron in top)
        cells = all_simplices(self.top)
        self.f_vector = tuple(len(layer) for layer in cells)
        edges, faces = cells[1], cells[2]
        self.edge_classes, edge_residual = classify_entities(
            self.positions, edges, edge_signature
        )
        self.face_classes, face_residual = classify_entities(
            self.positions, faces, face_signature
        )
        self.tetra_classes, tetra_residual = classify_entities(
            self.positions, self.top, tetra_signature
        )
        self.maximum_class_residual = max(
            edge_residual, face_residual, tetra_residual
        )

        self.edge_points = np.asarray([
            record["representative_points"] for record in self.edge_classes
        ])
        self.edge_weights = np.asarray([
            record["multiplicity"] for record in self.edge_classes
        ], dtype=float)
        self.edge_lengths = np.linalg.norm(
            self.edge_points[:, 0]-self.edge_points[:, 1], axis=1
        )
        self.face_points = np.asarray([
            record["representative_points"] for record in self.face_classes
        ])
        self.face_weights = np.asarray([
            record["multiplicity"] for record in self.face_classes
        ], dtype=float)
        self.face_areas = triangle_area(self.face_points)
        self.points = np.asarray([
            record["representative_points"] for record in self.tetra_classes
        ])
        self.tetra_weights = np.asarray([
            record["multiplicity"] for record in self.tetra_classes
        ], dtype=float)
        self.tetra_volumes = tetra_volume(self.points)
        self.volume_bar = float(np.dot(self.tetra_weights,
                                       self.tetra_volumes))
        self.s0 = float((2*np.pi**2/self.volume_bar)**(1/3))

        self.differences = self.points[:, :3, :]-self.points[:, 3:4, :]
        self.spatial_gram = np.einsum(
            "nai,nbi->nab", self.differences, self.differences
        )
        self.cross_shape = np.einsum(
            "nai,ni->na", self.differences, self.points[:, 3, :]
        )
        self.local_edge_lengths = np.stack([
            np.linalg.norm(self.points[:, left]-self.points[:, right], axis=1)
            for left, right in LOCAL_PAIRS
        ], axis=1)
        self.local_face_areas = np.stack([
            triangle_area(self.points[:, [j for j in range(4) if j != i], :])
            for i in range(4)
        ], axis=1)
        self.spatial_angles = self._spatial_angles()
        curvature_bar = (
            2*np.pi*np.dot(self.edge_weights, self.edge_lengths)
            - np.sum(self.tetra_weights[:, None]
                     * self.local_edge_lengths*self.spatial_angles)
        )
        self.curvature = float(self.s0*curvature_bar)
        self.mass = float(self.curvature/(8*np.pi))

    @staticmethod
    def _inner(left, inverse, right):
        return np.einsum("nki,nij,nkj->nk", left, inverse, right)

    def _spatial_angles(self):
        inverse = np.linalg.inv(self.spatial_gram)
        normals = np.zeros((len(self.points), 4, 3))
        normals[:, 0, 0] = 1
        normals[:, 1, 1] = 1
        normals[:, 2, 2] = 1
        normals[:, 3, :] = -1
        result = []
        for left, right in LOCAL_PAIRS:
            first, second = [value for value in range(4)
                             if value not in (left, right)]
            n1 = normals[:, first:first+1]
            n2 = normals[:, second:second+1]
            cosine = -self._inner(n1, inverse, n2)/np.sqrt(
                self._inner(n1, inverse, n1)
                * self._inner(n2, inverse, n2)
            )
            result.append(np.arccos(np.clip(cosine[:, 0], -1, 1)))
        return np.stack(result, axis=1)

    def _metric_data(self, s_minus, s_plus, rho):
        dtype = np.result_type(s_minus, s_plus, rho)
        metric = np.empty((len(self.points), 4, 4), dtype=dtype)
        metric[:, :3, :3] = s_minus*s_minus*self.spatial_gram
        cross = s_minus*(s_plus-s_minus)*self.cross_shape
        metric[:, :3, 3] = cross
        metric[:, 3, :3] = cross
        metric[:, 3, 3] = -rho
        inverse = np.linalg.inv(metric)
        normals = np.zeros((len(self.points), 4, 4), dtype=dtype)
        normals[:, 0, 0] = 1
        normals[:, 1, 1] = 1
        normals[:, 2, 2] = 1
        normals[:, 3, :3] = -1
        normals[:, 3, 3] = s_plus/s_minus-1
        bottom = np.zeros((len(self.points), 1, 4), dtype=dtype)
        bottom[:, :, 3] = 1
        return inverse, normals, bottom, -bottom

    def _angle(self, left, inverse, right, boundary=False):
        cross = self._inner(left, inverse, right)
        product = (self._inner(left, inverse, left)
                   * self._inner(right, inverse, right))
        denominator = (-1j*np.sqrt(-product.astype(complex)) if boundary
                       else np.sqrt(product.astype(complex)))
        return np.arccos(-cross/denominator)

    def gravitational_action(self, s_minus, s_plus, rho,
                             tetra_weight_delta=None):
        inverse, normals, bottom, top = self._metric_data(
            s_minus, s_plus, rho
        )
        weights = self.tetra_weights.copy()
        if tetra_weight_delta is not None:
            weights[tetra_weight_delta] += 1
        lateral = []
        for left, right in LOCAL_PAIRS:
            first, second = [value for value in range(4)
                             if value not in (left, right)]
            lateral.append(self._angle(
                normals[:, first:first+1], inverse,
                normals[:, second:second+1], False
            )[:, 0])
        lateral = np.stack(lateral, axis=1)
        bottom_angles = np.stack([
            self._angle(bottom, inverse, normals[:, index:index+1], True)[:, 0]
            for index in range(4)
        ], axis=1)
        top_angles = np.stack([
            self._angle(top, inverse, normals[:, index:index+1], True)[:, 0]
            for index in range(4)
        ], axis=1)

        delta = s_plus-s_minus
        edge_areas = (1j*(s_minus+s_plus)/2*self.edge_lengths
                      * np.sqrt(rho+delta*delta*self.edge_lengths**2/4))
        local_edge_areas = (
            1j*(s_minus+s_plus)/2*self.local_edge_lengths
            * np.sqrt(rho+delta*delta*self.local_edge_lengths**2/4)
        )
        global_face_area = np.dot(self.face_weights, self.face_areas)
        sum_hinges = (
            2*np.pi*np.dot(self.edge_weights, edge_areas)
            - np.sum(weights[:, None]*local_edge_areas*lateral)
            + np.pi*(s_minus*s_minus+s_plus*s_plus)*global_face_area
            - np.sum(weights[:, None]
                     * s_minus*s_minus*self.local_face_areas*bottom_angles)
            - np.sum(weights[:, None]
                     * s_plus*s_plus*self.local_face_areas*top_angles)
        )
        return -1j*sum_hinges

    def total_action(self, s_minus, s_plus, rho):
        return (self.gravitational_action(s_minus, s_plus, rho)
                - 8*np.pi*self.mass*np.sqrt(rho))

    def census(self):
        return {
            "f_vector": list(self.f_vector),
            "edge_classes": len(self.edge_classes),
            "face_classes": len(self.face_classes),
            "tetrahedron_classes": len(self.tetra_classes),
            "edge_multiplicity_sum": int(self.edge_weights.sum()),
            "face_multiplicity_sum": int(self.face_weights.sum()),
            "tetrahedron_multiplicity_sum": int(self.tetra_weights.sum()),
            "maximum_class_residual": self.maximum_class_residual,
            "volume_bar": self.volume_bar,
            "scale_for_unit_volume_radius": self.s0,
            "curvature": self.curvature,
            "mass": self.mass,
        }


def five_point_log_derivative(mesh, state, coordinate, step):
    values = [float(value) for value in state]

    def evaluate(offset):
        shifted = values.copy()
        shifted[coordinate] *= np.exp(offset*step)
        return mesh.total_action(*shifted)

    return (evaluate(-2)-8*evaluate(-1)+8*evaluate(1)-evaluate(2))/(12*step)


def richardson(values):
    values = np.asarray(values, dtype=float)
    first = (4*values[1:]-values[:-1])/3
    second = (16*first[1:]-first[:-1])/15
    return float(second[-1])


def affine_root(values):
    a0, a1 = SENTINELS[:2]
    beta = (values[a1]-values[a0])/(a1-a0)
    return -values[a0]/beta


def quadratic_root(values):
    a1, a2 = SENTINELS[1:3]
    matrix = np.array(((a1, a1*a1), (a2, a2*a2)))
    linear, quadratic = np.linalg.solve(
        matrix, np.array((values[a1], values[a2]))
    )
    root = -linear/quadratic
    a3 = SENTINELS[3]
    predicted = linear*a3+quadratic*a3*a3
    residual = abs(predicted-values[a3])/max(
        abs(predicted), abs(values[a3]), 1e-30
    )
    return float(root), float(residual)


def coefficient_route(mesh, seam_step, lapse_step):
    seam_tables = {a: [] for a in SENTINELS}
    for eta in ETAS_SEAM:
        rho = eta*eta
        previous = 0.5*five_point_log_derivative(
            mesh, (mesh.s0, mesh.s0, rho), 1, seam_step
        )
        for a in SENTINELS:
            s_plus = mesh.s0*np.exp(a*eta*eta)
            current = 0.5*five_point_log_derivative(
                mesh, (mesh.s0, s_plus, rho), 0, seam_step
            )
            seam_tables[a].append(float(np.real((previous+current)/eta)))
    seam_limits = {a: richardson(seam_tables[a]) for a in SENTINELS}
    seam_root = affine_root(seam_limits)

    lapse_tables = {a: [] for a in SENTINELS}
    for eta in ETAS_LAPSE:
        rho = eta*eta
        for a in SENTINELS:
            s_plus = mesh.s0*np.exp(a*eta*eta)
            derivative = five_point_log_derivative(
                mesh, (mesh.s0, s_plus, rho), 2, lapse_step
            )
            lapse_tables[a].append(float(np.real(derivative/eta**3)))
    lapse_limits = {a: richardson(lapse_tables[a]) for a in SENTINELS}
    lapse_root, quadratic_residual = quadratic_root(lapse_limits)
    static_residual = abs(lapse_limits[0.0])/max(
        max(abs(value) for value in lapse_limits.values()), 1e-30
    )
    return {
        "seam_coefficient": seam_root,
        "lapse_coefficient": lapse_root,
        "seam_lapse_difference": abs(seam_root-lapse_root),
        "lapse_quadratic_residual": quadratic_residual,
        "lapse_static_residual": static_residual,
        "seam_limits": {str(key): value for key, value in seam_limits.items()},
        "lapse_limits": {str(key): value for key, value in lapse_limits.items()},
    }


def coefficient_audit(mesh):
    primary = coefficient_route(mesh, SEAM_STEPS[0], LAPSE_STEPS[0])
    secondary = coefficient_route(mesh, SEAM_STEPS[1], LAPSE_STEPS[1])
    return {
        "coefficient": primary["seam_coefficient"],
        "primary": primary,
        "secondary": secondary,
        "seam_step_difference": abs(
            primary["seam_coefficient"]-secondary["seam_coefficient"]
        ),
        "lapse_step_difference": abs(
            primary["lapse_coefficient"]-secondary["lapse_coefficient"]
        ),
    }


def load_full_mesh_class():
    wanted = {"mesh_edges", "mesh_faces_with_counts", "triangle_area",
              "tetra_volume", "CellularMesh"}
    tree = ast.parse(ACTION_SOURCE.read_text(), filename=str(ACTION_SOURCE))
    selected = [node for node in tree.body
                if isinstance(node, (ast.FunctionDef, ast.ClassDef))
                and node.name in wanted]
    namespace = {
        "np": np, "Counter": Counter, "combinations": combinations,
        "hashlib": hashlib, "ETAS_G": (0.04, 0.02, 0.01, 0.005),
        "ETAS_F": (0.04, 0.02, 0.01),
        "A_SENTINELS": (0.0, -1.0, -2.0, -3.0),
        "PRIMARY_DERIVATIVE_STEP": 2e-5,
        "SECONDARY_DERIVATIVE_STEP": 1e-5,
        "PRIMARY_LAPSE_DERIVATIVE_STEP": 4e-3,
        "SECONDARY_LAPSE_DERIVATIVE_STEP": 8e-3,
    }
    exec(compile(ast.Module(body=selected, type_ignores=[]),
                 str(ACTION_SOURCE), "exec"), namespace)
    return namespace["CellularMesh"]


print("="*78)
print("ADVERSARIAL ORBIT-COMPRESSED CANONICAL ACCELERATION")
print("="*78)

actual_hashes = {path.name: digest(path) for path in INPUT_HASHES}
expected_hashes = {path.name: value for path, value in INPUT_HASHES.items()}
provenance_ok = check(
    "the audit has exact frozen primary provenance",
    actual_hashes == expected_hashes and PROTOCOL_COMMIT == "3b7bd6c"
    and PROTOCOL_CORRECTION_COMMIT == "1c3318d",
    str(actual_hashes),
)
primary = json.loads(PRIMARY_BLIND.read_text())

vertices, adjacency, _ = build_600cell()
vertices = vertices/np.linalg.norm(vertices, axis=1)[:, None]
coarse_top = networkx_tetrahedra(adjacency)
check(
    "NetworkX independently reconstructs the 600-cell boundary",
    len(coarse_top) == 600 and {len(tet) for tet in coarse_top} == {4},
)
coarse_cells = all_simplices(coarse_top)
base_positions, base_top = projected_barycentric(
    vertices, coarse_cells, coarse_top
)
fine_positions, fine_top = direct_rank_split(base_positions, base_top)

carriers = {
    "regular_600cell_control": (vertices, coarse_top),
    "projected_barycentric": (base_positions, base_top),
    "projected_rank_edgewise_2": (fine_positions, fine_top),
}
reduced = {}
for name, (positions, top) in carriers.items():
    print(f"[INFO] classifying {name}", flush=True)
    reduced[name] = ReducedMesh(name, positions, top)
    census = reduced[name].census()
    check(
        f"{name} has an exact intrinsic class census",
        census["maximum_class_residual"] < 2e-12
        and census["edge_multiplicity_sum"] == census["f_vector"][1]
        and census["face_multiplicity_sum"] == census["f_vector"][2]
        and census["tetrahedron_multiplicity_sum"] == census["f_vector"][3],
        str(census),
    )

FullMesh = load_full_mesh_class()
compression_errors = {}
for name in ("projected_barycentric", "projected_rank_edgewise_2"):
    positions, top = carriers[name]
    full = FullMesh(name+"_heldout", positions, top)
    mesh = reduced[name]
    states = (
        (mesh.s0, mesh.s0, 0.02**2),
        (mesh.s0, mesh.s0*0.97, 0.02),
        (mesh.s0*1.01, mesh.s0*0.985, 0.015),
    )
    errors = []
    for state_index, state in enumerate(states):
        if state_index == 0:
            expected = full.gravitational_action(*state)[0]
            observed = mesh.gravitational_action(*state)
        else:
            expected = full.total_action(*state)
            observed = mesh.total_action(*state)
        errors.append(float(abs(observed-expected)/max(
            abs(observed), abs(expected), 1e-30
        )))
    compression_errors[name] = errors
    check(
        f"{name} compression matches complete held-out actions",
        max(errors) < 2e-8,
        f"relative errors={errors}",
    )

fine = reduced["projected_rank_edgewise_2"]
largest_class = int(np.argmax(fine.tetra_weights))
heldout = (fine.s0, fine.s0*0.97, 0.02)
correct = fine.gravitational_action(*heldout)
wrong = fine.gravitational_action(*heldout,
                                  tetra_weight_delta=largest_class)
negative_difference = float(abs(wrong-correct)/max(abs(correct), 1e-30))
check(
    "the wrong-multiplicity negative control is detected",
    negative_difference > 1e-8,
    f"relative difference={negative_difference:.3e}",
)

audits = {}
for name, mesh in reduced.items():
    print(f"[INFO] extracting coefficient for {name}", flush=True)
    audits[name] = coefficient_audit(mesh)

regular_error = abs(
    audits["regular_600cell_control"]["coefficient"]-REGULAR_CONTROL
)
check(
    "the independent derivative route reproduces the regular control",
    regular_error < 2e-5,
    f"a={audits['regular_600cell_control']['coefficient']:.12g}, "
    f"error={regular_error:.3e}",
)

coefficient_errors = {}
route_ok = True
for name in ("projected_barycentric", "projected_rank_edgewise_2"):
    audit = audits[name]
    target = primary["blind_coefficients"][name]
    coefficient_errors[name] = abs(audit["coefficient"]-target)
    route_ok &= bool(
        audit["primary"]["seam_lapse_difference"] < 2e-5
        and audit["secondary"]["seam_lapse_difference"] < 2e-5
        and audit["seam_step_difference"] < 2e-5
        and audit["lapse_step_difference"] < 2e-5
        and coefficient_errors[name] < 2e-5
        and audit["primary"]["lapse_quadratic_residual"] < 2e-5
        and audit["primary"]["lapse_static_residual"] < 2e-5
    )
check(
    "both orbit-compressed coefficients agree by every frozen route",
    route_ok,
    f"primary differences={coefficient_errors}",
)

all_ok = bool(
    provenance_ok and route_ok and regular_error < 2e-5
    and negative_difference > 1e-8
    and max(max(values) for values in compression_errors.values()) < 2e-8
    and all(mesh.maximum_class_residual < 2e-12
            for mesh in reduced.values())
)
outcome = (
    "ADVERSARIAL_CANONICAL_ACCELERATION_CORROBORATED" if all_ok else
    "ADVERSARIAL_CANONICAL_ACCELERATION_DISAGREEMENT"
)
check(
    "the adversarial hierarchy assigns exactly one outcome",
    outcome in {
        "ADVERSARIAL_CANONICAL_ACCELERATION_CORROBORATED",
        "ADVERSARIAL_CANONICAL_ACCELERATION_DISAGREEMENT",
    },
    outcome,
)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "protocol_correction_commit": PROTOCOL_CORRECTION_COMMIT,
    "input_sha256": actual_hashes,
    "method": "intrinsic-shape census plus real five-point derivatives",
    "censuses": {name: mesh.census() for name, mesh in reduced.items()},
    "compression_relative_errors": compression_errors,
    "negative_control_relative_difference": negative_difference,
    "regular_control_target": REGULAR_CONTROL,
    "regular_control_error": regular_error,
    "coefficient_audits": audits,
    "primary_coefficient_differences": coefficient_errors,
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print("="*78)
print(f"Primary coefficient differences: {coefficient_errors}")
print(f"Outcome: {outcome}")
print(f"{passed}/{tests} checks passed")
sys.exit(0 if passed == tests and all_ok else 1)

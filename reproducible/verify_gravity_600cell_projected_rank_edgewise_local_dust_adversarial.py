#!/usr/bin/env python3
"""Consistent-mass-matrix audit of conditional local P1 dust weights."""

from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import networkx as nx
import numpy as np
from scipy import sparse
import sympy as sy


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons import build_600cell  # noqa: E402


PRIMARY = HERE / "gravity_600cell_projected_rank_edgewise_local_dust.json"
OUTPUT = HERE / "gravity_600cell_projected_rank_edgewise_local_dust_adversarial.json"
PRIMARY_SHA256 = "53463e5271301ae41eb26564875d26991ddea8024a9e09ae3c302d428ad39779"
PROTOCOL_COMMIT = "0257417"
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


def all_simplices(top):
    return tuple(
        tuple(sorted({tuple(sorted(face)) for tetrahedron in top
                      for face in combinations(tetrahedron, degree+1)}))
        for degree in range(4)
    )


def networkx_tetrahedra(adjacency):
    graph = nx.from_numpy_array(adjacency)
    return tuple(sorted(tuple(sorted(clique))
                        for clique in nx.find_cliques(graph)))


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
    keys = set()
    for chamber in base_top:
        names = {f"v{rank}": (chamber[rank], chamber[rank])
                 for rank in range(4)}
        for left, right in combinations(range(4), 2):
            names[f"m{left}{right}"] = tuple(sorted(
                (chamber[left], chamber[right])
            ))
        for child in DIRECT_CHILDREN:
            tetrahedron = tuple(sorted(names[name] for name in child))
            top_keys.append(tetrahedron)
            keys.update(tetrahedron)
    ordered = tuple(sorted(keys))
    index = {key: position for position, key in enumerate(ordered)}
    raw = np.asarray([base_positions[list(key)].mean(axis=0)
                      for key in ordered])
    positions = raw/np.linalg.norm(raw, axis=1)[:, None]
    top = tuple(tuple(sorted(index[key] for key in tetrahedron))
                for tetrahedron in top_keys)
    return positions, top


def gram_volumes(positions, top):
    tetrahedra = np.asarray(top, dtype=np.int32)
    values = []
    for start in range(0, len(tetrahedra), 20000):
        points = positions[tetrahedra[start:start+20000]]
        edges = points[:, 1:, :]-points[:, :1, :]
        gram = np.einsum("nai,nbi->nab", edges, edges)
        values.append(np.sqrt(np.maximum(np.linalg.det(gram), 0.0))/6)
    return np.concatenate(values)


def cayley_menger_volumes(positions, top):
    tetrahedra = np.asarray(top, dtype=np.int32)
    values = []
    for start in range(0, len(tetrahedra), 12000):
        points = positions[tetrahedra[start:start+12000]]
        differences = points[:, :, None, :]-points[:, None, :, :]
        squared = np.einsum("nija,nija->nij", differences, differences)
        matrix = np.ones((len(points), 5, 5), dtype=float)
        matrix[:, 0, 0] = 0
        matrix[:, 1:, 1:] = squared
        values.append(np.sqrt(np.maximum(np.linalg.det(matrix)/288, 0.0)))
    return np.concatenate(values)


def direct_weights(volumes, top, vertex_count):
    tetrahedra = np.asarray(top, dtype=np.int32)
    return np.bincount(
        tetrahedra.reshape(-1), weights=np.repeat(volumes/4, 4),
        minlength=vertex_count,
    )


def consistent_matrix(volumes, top, vertex_count):
    tetrahedra = np.asarray(top, dtype=np.int32)
    rows = np.repeat(tetrahedra, 4, axis=1).reshape(-1)
    cols = np.tile(tetrahedra, (1, 4)).reshape(-1)
    template = np.ones((4, 4), dtype=float)+np.eye(4)
    data = (volumes[:, None, None]*template[None, :, :]/20).reshape(-1)
    matrix = sparse.coo_matrix(
        (data, (rows, cols)), shape=(vertex_count, vertex_count)
    ).tocsr()
    return matrix


def weight_digest(weights):
    return sha256(np.asarray(weights, dtype="<f8").tobytes()).hexdigest()


print("="*78)
print("ADVERSARIAL CONSISTENT-MATRIX LOCAL P1 DUST")
print("="*78)

primary_hash = sha256(PRIMARY.read_bytes()).hexdigest()
primary = json.loads(PRIMARY.read_text())
provenance_ok = check(
    "the audit has exact frozen primary provenance",
    primary_hash == PRIMARY_SHA256 and PROTOCOL_COMMIT == "0257417"
    and primary["outcome"] == "P1_LOCAL_DUST_WEIGHTS_DERIVED_CONDITIONALLY",
    f"primary={primary_hash}",
)

volume_symbol = sy.symbols("V", positive=True)
local_matrix = volume_symbol/sy.Integer(20)*(
    sy.ones(4, 4)+sy.eye(4)
)
symbolic_rows = tuple(sy.simplify(sum(local_matrix.row(index)))
                      for index in range(4))
check(
    "the consistent local P1 matrix has exact row sum V/4",
    symbolic_rows == (volume_symbol/4,)*4,
    str(symbolic_rows),
)

regular = np.asarray(((1, 1, 1, 0), (1, -1, -1, 0),
                      (-1, 1, -1, 0), (-1, -1, 1, 0)), dtype=float)
regular_top = ((0, 1, 2, 3),)
regular_gram = gram_volumes(regular, regular_top)[0]
regular_cayley = cayley_menger_volumes(regular, regular_top)[0]
regular_error = abs(regular_gram-regular_cayley)/regular_gram
check(
    "Cayley--Menger reproduces the regular-tetrahedron Gram volume",
    regular_error < 2e-14,
    f"Gram={regular_gram:.15g}, Cayley={regular_cayley:.15g}, "
    f"relative={regular_error:.3e}",
)

vertices, adjacency, _ = build_600cell()
vertices = vertices/np.linalg.norm(vertices, axis=1)[:, None]
coarse_top = networkx_tetrahedra(adjacency)
check(
    "NetworkX independently finds the 600 source tetrahedra",
    len(coarse_top) == 600 and {len(tet) for tet in coarse_top} == {4},
)
coarse_cells = all_simplices(coarse_top)
base_positions, base_top = projected_barycentric(
    vertices, coarse_cells, coarse_top
)
fine_positions, fine_top = direct_rank_split(base_positions, base_top)

levels = {
    "regular_600cell_control": (vertices, coarse_top),
    "projected_barycentric": (base_positions, base_top),
    "projected_rank_edgewise_2": (fine_positions, fine_top),
}
records = {}
all_levels_ok = True
for name, (positions, top) in levels.items():
    print(f"[INFO] assembling {name}", flush=True)
    gram = gram_volumes(positions, top)
    cayley = cayley_menger_volumes(positions, top)
    gram_weights = direct_weights(gram, top, len(positions))
    matrix = consistent_matrix(cayley, top, len(positions))
    row_weights = np.asarray(matrix.sum(axis=1)).ravel()
    difference = matrix-matrix.T
    symmetry_residual = (0.0 if difference.nnz == 0 else
                         float(np.max(np.abs(difference.data))))
    pointwise_relative = float(np.max(
        np.abs(row_weights-gram_weights)
        / np.maximum(np.abs(gram_weights), 1e-300)
    ))
    minimum_local_eigenvalue = float(cayley.min()/20)

    if name == "regular_600cell_control":
        uniform_residual = float((row_weights.max()-row_weights.min())
                                 / row_weights.mean())
        control_ok = check(
            "the regular 600-cell row-sum masses are uniform",
            uniform_residual < 2e-10,
            f"relative spread={uniform_residual:.3e}",
        )
        frozen_digest_ok = True
        frozen_scalar_ok = True
        negative_discrepancy = None
    else:
        frozen = primary["levels"][name]
        frozen_digest_ok = weight_digest(gram_weights) == frozen["weight_sha256"]
        scalar_errors = (
            abs(gram_weights.sum()-frozen["weight_sum"])
            / frozen["weight_sum"],
            abs(gram_weights.min()-frozen["weight_minimum"])
            / frozen["weight_minimum"],
            abs(gram_weights.max()-frozen["weight_maximum"])
            / frozen["weight_maximum"],
        )
        frozen_scalar_ok = max(scalar_errors) < 2e-10
        uniform = gram_weights.sum()/len(gram_weights)
        negative_discrepancy = float(np.max(
            np.abs(uniform-gram_weights)/gram_weights
        ))
        control_ok = check(
            f"{name} rejects globally uniform masses as affine-exact",
            negative_discrepancy > 0.1,
            f"maximum P1 hat-function discrepancy={negative_discrepancy:.6g}",
        )

    assembly_ok = check(
        f"{name} consistent row sums reproduce direct Gram weights",
        pointwise_relative < 2e-10 and symmetry_residual < 2e-14
        and minimum_local_eigenvalue > 0 and row_weights.min() > 0
        and frozen_digest_ok and frozen_scalar_ok,
        f"pointwise={pointwise_relative:.3e}, symmetry={symmetry_residual:.3e}, "
        f"digest={frozen_digest_ok}",
    )
    all_levels_ok &= control_ok and assembly_ok
    records[name] = {
        "vertices": len(positions),
        "tetrahedra": len(top),
        "gram_total_volume": float(gram.sum()),
        "cayley_total_volume": float(cayley.sum()),
        "consistent_matrix_nnz": int(matrix.nnz),
        "matrix_symmetry_residual": symmetry_residual,
        "minimum_local_matrix_eigenvalue": minimum_local_eigenvalue,
        "maximum_pointwise_relative_difference": pointwise_relative,
        "gram_weight_sha256": weight_digest(gram_weights),
        "row_sum_weight_sha256": weight_digest(row_weights),
        "uniform_negative_control_discrepancy": negative_discrepancy,
    }

all_ok = bool(provenance_ok and all_levels_ok and regular_error < 2e-14)
outcome = (
    "ADVERSARIAL_P1_LOCAL_DUST_CORROBORATED" if all_ok else
    "ADVERSARIAL_P1_LOCAL_DUST_DISAGREEMENT"
)
check(
    "the adversarial hierarchy assigns exactly one outcome",
    outcome in {"ADVERSARIAL_P1_LOCAL_DUST_CORROBORATED",
                "ADVERSARIAL_P1_LOCAL_DUST_DISAGREEMENT"},
    outcome,
)

artifact = {
    "protocol_commit": PROTOCOL_COMMIT,
    "primary_sha256": primary_hash,
    "method": "Cayley-Menger volumes plus consistent P1 mass-matrix row sums",
    "symbolic_local_row_sums": [str(value) for value in symbolic_rows],
    "regular_tetrahedron_volume_relative_error": regular_error,
    "levels": records,
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print("="*78)
print(f"Outcome: {outcome}")
print(f"{passed}/{tests} checks passed")
sys.exit(0 if passed == tests and all_ok else 1)

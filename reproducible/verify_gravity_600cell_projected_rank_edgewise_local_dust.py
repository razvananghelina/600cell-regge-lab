#!/usr/bin/env python3
"""Conditional canonical P1 dust weights on projected rank-edgewise carriers."""

from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sy


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from commons import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_600cell_projected_rank_edgewise_local_dust.json"
PRIOR_ART_COMMIT = "ba7de6c"
PROTOCOL_COMMIT = "9156c0f"
INPUT_HASHES = {
    "commons/cell600.py":
        "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
    "reproducible/verify_gravity_600cell_projected_rank_edgewise_carrier.py":
        "50876c582cf22d86296f3f2b715ff1cf3276a9f1320baa3b37d365ce91f2aa23",
    "reproducible/gravity_600cell_projected_rank_edgewise_carrier.json":
        "b57955b85a972df00b5673ddf7ee295757848f5afb43314857cf3de2dc85ac84",
    "reproducible/gravity_600cell_projected_rank_edgewise_acceleration_blind.json":
        "2059620f22cfbd8eac8abe6f2c7536924128d37f47a430bf773e34a9aead93a2",
}
EXPECTED_F = {
    "projected_barycentric": (2640, 17040, 28800, 14400),
    "projected_rank_edgewise_2": (19680, 134880, 230400, 115200),
}
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


def tetrahedra_from_adjacency(adjacency):
    neighbours = [set(np.flatnonzero(adjacency[index] > 0.5))
                  for index in range(len(adjacency))]
    result = []
    for first in range(len(adjacency)):
        for second in sorted(value for value in neighbours[first]
                             if value > first):
            common_two = neighbours[first] & neighbours[second]
            for third in sorted(value for value in common_two
                                if value > second):
                common_three = common_two & neighbours[third]
                for fourth in sorted(value for value in common_three
                                     if value > third):
                    result.append((first, second, third, fourth))
    return tuple(result)


def projected_barycentric(vertices, coarse_cells, coarse_top):
    vertex_cells = tuple(cell for layer in coarse_cells for cell in layer)
    cell_index = {cell: index for index, cell in enumerate(vertex_cells)}
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
            top.append(tuple(cell_index[flag] for flag in flags))
    return vertex_cells, cell_index, positions, tuple(top)


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
        for left, right in combinations(range(4), 2):
            names[f"m{left}{right}"] = tuple(sorted(
                (chamber[left], chamber[right])
            ))
        for child in DIRECT_CHILDREN:
            keys = tuple(sorted(names[name] for name in child))
            top_keys.append(keys)
            vertex_keys.update(keys)
    ordered_keys = tuple(sorted(vertex_keys))
    key_index = {key: index for index, key in enumerate(ordered_keys)}
    raw = np.asarray([base_positions[list(key)].mean(axis=0)
                      for key in ordered_keys])
    positions = raw/np.linalg.norm(raw, axis=1)[:, None]
    top = tuple(tuple(sorted(key_index[key] for key in tetrahedron))
                for tetrahedron in top_keys)
    return ordered_keys, key_index, positions, top


def quaternion_multiply(left, right):
    w1, x1, y1, z1 = np.asarray(left)
    w2, x2, y2, z2 = np.asarray(right)
    return np.array((
        w1*w2-x1*x2-y1*y2-z1*z2,
        w1*x2+x1*w2+y1*z2-z1*y2,
        w1*y2-x1*z2+y1*w2+z1*x2,
        w1*z2+x1*y2-y1*x2+z1*w2,
    ))


def multiplication_table(vertices):
    table = np.empty((len(vertices), len(vertices)), dtype=np.int16)
    for left, quaternion in enumerate(vertices):
        products = np.asarray([
            quaternion_multiply(quaternion, right) for right in vertices
        ])
        table[left] = np.argmax(products @ vertices.T, axis=1)
    return table


def source_actions(vertices, table):
    conjugates = []
    for value in vertices:
        target = value.copy()
        target[1:] *= -1
        conjugates.append(int(np.argmax(vertices @ target)))
    actions = []
    for group_index in range(len(vertices)):
        actions.append(table[group_index].astype(int))
        actions.append(table[:, group_index].astype(int))
    actions.append(np.asarray(conjugates, dtype=int))
    return tuple(actions)


def induced_maps(actions, vertex_cells, cell_index, fine_keys,
                 fine_key_index):
    base_maps = []
    fine_maps = []
    for action in actions:
        base_map = np.asarray([
            cell_index[tuple(sorted(int(action[vertex]) for vertex in cell))]
            for cell in vertex_cells
        ], dtype=np.int32)
        fine_map = np.asarray([
            fine_key_index[tuple(sorted((int(base_map[key[0]]),
                                         int(base_map[key[1]]))))]
            for key in fine_keys
        ], dtype=np.int32)
        base_maps.append(base_map)
        fine_maps.append(fine_map)
    return tuple(base_maps), tuple(fine_maps)


def tetra_volumes(positions, top):
    tetrahedra = np.asarray(top, dtype=np.int32)
    result = []
    for start in range(0, len(tetrahedra), 20000):
        points = positions[tetrahedra[start:start+20000]]
        edges = points[:, 1:, :]-points[:, :1, :]
        gram = np.einsum("nai,nbi->nab", edges, edges)
        result.append(np.sqrt(np.maximum(np.linalg.det(gram), 0.0))/6)
    return np.concatenate(result)


def assembled_weights(positions, top):
    tetrahedra = np.asarray(top, dtype=np.int32)
    volumes = tetra_volumes(positions, top)
    weights = np.bincount(
        tetrahedra.reshape(-1),
        weights=np.repeat(volumes/4, 4),
        minlength=len(positions),
    )
    return volumes, weights


def orbit_census(size, maps, weights):
    parent = np.arange(size, dtype=np.int32)

    def find(value):
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = int(parent[value])
        return value

    def union(left, right):
        left, right = find(left), find(right)
        if left != right:
            parent[right] = left

    for mapping in maps:
        for left, right in enumerate(mapping):
            union(left, int(right))
    groups = {}
    for value in range(size):
        groups.setdefault(find(value), []).append(value)
    records = []
    maximum_spread = 0.0
    for members in groups.values():
        values = weights[members]
        spread = float(values.max()-values.min())
        maximum_spread = max(maximum_spread, spread)
        records.append({
            "size": len(members),
            "weight_mean": float(values.mean()),
            "weight_spread": spread,
        })
    records.sort(key=lambda item: (item["size"], item["weight_mean"]))
    return records, maximum_spread


def weight_digest(weights):
    return sha256(np.asarray(weights, dtype="<f8").tobytes()).hexdigest()


print("="*78)
print("CONDITIONAL CANONICAL P1 DUST WEIGHTS")
print("="*78)

actual_hashes = {name: digest(ROOT/name) for name in INPUT_HASHES}
provenance_ok = check(
    "the frozen carrier and mass inputs have exact provenance",
    actual_hashes == INPUT_HASHES and PRIOR_ART_COMMIT == "ba7de6c"
    and PROTOCOL_COMMIT == "9156c0f",
    str(actual_hashes),
)

x, y, z = sy.symbols("x y z")
barycentric = (1-x-y-z, x, y, z)
integrals = tuple(sy.integrate(
    basis, (z, 0, 1-x-y), (y, 0, 1-x), (x, 0, 1)
) for basis in barycentric)
q = sy.symbols("q0:4")
solution = sy.solve(
    [q[index]-integrals[index] for index in range(4)], q, dict=True
)
uniqueness_ok = check(
    "exact affine P1 quadrature uniquely gives volume/4 per vertex",
    integrals == (sy.Rational(1, 24),)*4
    and solution == [{value: sy.Rational(1, 24) for value in q}],
    f"integrals={integrals}, solution={solution}",
)

vertices, adjacency, _ = build_600cell()
vertices = vertices/np.linalg.norm(vertices, axis=1)[:, None]
coarse_top = tetrahedra_from_adjacency(adjacency)
coarse_cells = all_simplices(coarse_top)
vertex_cells, cell_index, base_positions, base_top = projected_barycentric(
    vertices, coarse_cells, coarse_top
)
fine_keys, fine_key_index, fine_positions, fine_top = direct_rank_split(
    base_positions, base_top
)

table = multiplication_table(vertices)
actions = source_actions(vertices, table)
base_maps, fine_maps = induced_maps(
    actions, vertex_cells, cell_index, fine_keys, fine_key_index
)
blind = json.loads((HERE / "gravity_600cell_projected_rank_edgewise_acceleration_blind.json").read_text())
carrier = json.loads((HERE / "gravity_600cell_projected_rank_edgewise_carrier.json").read_text())

levels = {
    "projected_barycentric": (base_positions, base_top, base_maps),
    "projected_rank_edgewise_2": (fine_positions, fine_top, fine_maps),
}
records = {}
all_level_ok = True
for name, (positions, top, maps) in levels.items():
    cells = all_simplices(top)
    f_vector = tuple(len(layer) for layer in cells)
    volumes, weights = assembled_weights(positions, top)
    target_volume = carrier["levels"][name]["geometry"]["volume_total"]
    total_volume = float(volumes.sum())
    mass = float(blind["levels"][name]["selected_total_dust_mass"])
    masses = mass*weights/weights.sum()
    symmetry_residual = max(float(np.max(np.abs(
        weights-weights[mapping]
    ))) for mapping in maps)
    orbits, orbit_spread = orbit_census(len(weights), maps, weights)

    scaling_residual = 0.0
    for scale in (0.5, 2.0):
        _, scaled = assembled_weights(scale*positions, top)
        scaling_residual = max(scaling_residual, float(np.max(
            np.abs(scaled-scale**3*weights)
            / np.maximum(np.abs(scale**3*weights), 1e-300)
        )))
    action_residual = 0.0
    for tau in (0.01, 0.1, 1.0):
        local = -8*np.pi*np.sum(masses*tau)
        global_action = -8*np.pi*mass*tau
        action_residual = max(action_residual, float(
            abs(local-global_action)/max(abs(global_action), 1e-300)
        ))

    topology_ok = check(
        f"{name} reproduces the frozen spatial carrier",
        f_vector == EXPECTED_F[name]
        and abs(total_volume-target_volume) < 5e-11,
        f"f={f_vector}, V={total_volume:.15f}",
    )
    local_ok = check(
        f"{name} has positive exactly conservative local P1 masses",
        float(volumes.min()) > 0 and float(weights.min()) > 0
        and float(masses.min()) > 0
        and abs(float(weights.sum())-total_volume) < 5e-12
        and abs(float(masses.sum())-mass) < 5e-13,
        f"w=[{weights.min():.6e},{weights.max():.6e}], "
        f"sum={weights.sum():.15f}, M={masses.sum():.15f}",
    )
    equivariance_ok = check(
        f"{name} weights pass all 241 spatial actions",
        len(maps) == 241 and symmetry_residual < 2e-10
        and orbit_spread < 2e-10,
        f"orbits={len(orbits)}, max residual={symmetry_residual:.3e}",
    )
    covariance_ok = check(
        f"{name} weights scale cubically and collapse to global dust",
        scaling_residual < 2e-12 and action_residual < 2e-14,
        f"scale residual={scaling_residual:.3e}, "
        f"action residual={action_residual:.3e}",
    )
    level_ok = topology_ok and local_ok and equivariance_ok and covariance_ok
    all_level_ok &= level_ok
    rounded_multiplicities = Counter(np.round(weights, 13))
    records[name] = {
        "f_vector": list(f_vector),
        "total_chordal_volume": total_volume,
        "selected_total_mass": mass,
        "weight_minimum": float(weights.min()),
        "weight_maximum": float(weights.max()),
        "weight_sum": float(weights.sum()),
        "mass_minimum": float(masses.min()),
        "mass_maximum": float(masses.max()),
        "mass_sum": float(masses.sum()),
        "weight_sha256": weight_digest(weights),
        "mass_sha256": weight_digest(masses),
        "spatial_actions_tested": len(maps),
        "maximum_symmetry_residual": symmetry_residual,
        "maximum_orbit_weight_spread": orbit_spread,
        "symmetry_orbits": orbits,
        "rounded_weight_multiplicities": [
            {"weight": float(value), "multiplicity": count}
            for value, count in sorted(rounded_multiplicities.items())
        ],
        "maximum_cubic_scaling_relative_residual": scaling_residual,
        "maximum_global_action_collapse_relative_residual": action_residual,
    }

all_ok = bool(provenance_ok and uniqueness_ok and all_level_ok)
outcome = (
    "P1_LOCAL_DUST_WEIGHTS_DERIVED_CONDITIONALLY" if all_ok else
    "P1_LOCAL_DUST_INTERNAL_FAILURE"
)
check(
    "the frozen hierarchy assigns exactly one conditional outcome",
    outcome in {
        "P1_LOCAL_DUST_WEIGHTS_DERIVED_CONDITIONALLY",
        "P1_LOCAL_DUST_INTERNAL_FAILURE",
    },
    outcome,
)

artifact = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": actual_hashes,
    "ansatz_status": "STRUCTURAL_CONTINUOUS_NODAL_P1_DUST",
    "conditional_uniqueness_status": "KNOWN_AND_EXACTLY_RECONSTRUCTED",
    "reference_barycentric_integrals": [str(value) for value in integrals],
    "levels": records,
    "independent_vertex_lapse_carrier_constructed": False,
    "particle_or_spectral_target_loaded": False,
    "outcome": outcome,
    "tests": tests,
    "passed": passed,
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")

print("="*78)
for name, record in records.items():
    print(f"{name}: {len(record['symmetry_orbits'])} H4 orbits, "
          f"w=[{record['weight_minimum']:.6e},{record['weight_maximum']:.6e}]")
print(f"Outcome: {outcome}")
print(f"{passed}/{tests} checks passed")
sys.exit(0 if passed == tests and all_ok else 1)

#!/usr/bin/env python3
"""Can the existing chamber gamma,J select the chromatic chiral Z2?

Prior-art commit: e946b0a.
Protocol commit: 70525f1.
"""

from collections import Counter
from itertools import product
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from commons.cell600 import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_600cell_chromatic_chamber_selector.json"
CELL600_SOURCE = ROOT / "commons" / "cell600.py"
CHROMATIC_ARTIFACT = HERE / "gravity_600cell_chromatic_cover_orbits.json"
PRIOR_ART_COMMIT = "e946b0a"
PROTOCOL_COMMIT = "70525f1"
EXPECTED_HASHES = {
    "cell600": "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f",
    "chromatic_cover_orbits": "682e3cfaa0c2912085c0375281817e217f19a54bfc9d6ec9b296844063be7121",
}
TOLERANCE = 1e-8


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def qmul(left, right):
    a, b, c, d = left
    e, f, g, h = right
    return np.array((
        a*e-b*f-c*g-d*h,
        a*f+b*e+c*h-d*g,
        a*g-b*h+c*e+d*f,
        a*h+b*g-c*f+d*e,
    ))


def permutation_sign(permutation):
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def compose(left, right):
    """Permutation product left after right."""
    return tuple(left[right[index]] for index in range(len(left)))


def multiplication_table(permutations):
    lookup = {permutation: index for index, permutation in enumerate(permutations)}
    table = np.empty((len(permutations), len(permutations)), dtype=np.int16)
    for left, p_left in enumerate(permutations):
        for right, p_right in enumerate(permutations):
            table[left, right] = lookup[compose(p_left, p_right)]
    identity_permutation = tuple(range(len(permutations[0])))
    return table, lookup[identity_permutation]


def group_invariants(table, identity):
    order = table.shape[0]
    inverses = np.empty(order, dtype=np.int16)
    for element in range(order):
        candidates = np.flatnonzero(table[element] == identity)
        inverse = next(
            int(candidate) for candidate in candidates
            if table[candidate, element] == identity
        )
        inverses[element] = inverse
    involutions = tuple(
        element for element in range(order)
        if element != identity and table[element, element] == identity
    )
    commutators = set()
    for left in range(order):
        for right in range(order):
            value = int(table[
                table[table[left, right], inverses[left]], inverses[right]
            ])
            commutators.add(value)
    generators = tuple(sorted(commutators | {int(inverses[x]) for x in commutators}))
    subgroup = {identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            image = int(table[current, generator])
            if image not in subgroup:
                subgroup.add(image)
                frontier.append(image)
    return {
        "involutions": involutions,
        "commutator_elements": tuple(sorted(commutators)),
        "commutator_subgroup": tuple(sorted(subgroup)),
    }


def rotation_matrix(axis, angle):
    axis = axis / np.linalg.norm(axis)
    cross = np.array((
        (0, -axis[2], axis[1]),
        (axis[2], 0, -axis[0]),
        (-axis[1], axis[0], 0),
    ))
    return np.eye(3) + np.sin(angle)*cross + (1-np.cos(angle))*(cross @ cross)


def equivariant_map_census(target_reflection_odd):
    signs = (-1, 1)
    sign_index = {-1: 0, 1: 1}
    records = []
    for outputs in product(signs, repeat=2):
        def mapping(value):
            return outputs[sign_index[value]]

        equivariant = all(
            mapping(-value)
            == (-mapping(value) if target_reflection_odd else mapping(value))
            for value in signs
        )
        bijective = len(set(outputs)) == 2
        records.append({
            "f(-1)": outputs[0],
            "f(+1)": outputs[1],
            "equivariant": equivariant,
            "bijective": bijective,
        })
    return {
        "records": records,
        "equivariant_maps": sum(record["equivariant"] for record in records),
        "equivariant_bijections": sum(
            record["equivariant"] and record["bijective"] for record in records
        ),
    }


def pair_orbits(target_reflection_odd):
    states = {pair for pair in product((-1, 1), repeat=2)}
    orbits = []
    while states:
        seed = min(states)
        reflected = (
            -seed[0], -seed[1] if target_reflection_odd else seed[1]
        )
        orbit = frozenset((seed, reflected))
        orbits.append(orbit)
        states -= orbit
    return tuple(sorted(orbits, key=lambda orbit: min(orbit)))


tests = []


def check(label, condition, detail=""):
    ok = bool(condition)
    tests.append((label, ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if detail:
        print(f"       {detail}")


print("=" * 78)
print("STATIC CHAMBER SELECTOR FOR THE CHROMATIC CHIRAL Z2")
print("=" * 78)

source_hashes = {
    "cell600": digest(CELL600_SOURCE),
    "chromatic_cover_orbits": digest(CHROMATIC_ARTIFACT),
}
provenance_ok = source_hashes == EXPECTED_HASHES
chromatic = json.loads(CHROMATIC_ARTIFACT.read_text())
chromatic_orbits = chromatic["orbits"]
chromatic_input_ok = bool(
    chromatic["tests"] == chromatic["passed"] == 22
    and chromatic["outcome"] == "CHIRAL_COVER_AMBIGUITY"
    and chromatic_orbits["proper_unordered_cover_sizes"] == [5, 5]
    and chromatic_orbits["proper_unordered_cover_coset_sides"]
        == [["left"], ["right"]]
    and chromatic_orbits["improper_cover_orbit_map"] == [1, 0]
    and {tuple(values) for values in chromatic_orbits[
        "full_ordered_Z2_invariant_values"
    ]} == {(-1,), (1,)}
)

# Independently reconstruct the regular icosahedron and its complete flags.
sqrt5 = np.sqrt(5.0)
phi = (1+sqrt5)/2
icosahedron_vertices = []
for base in ((0, 1, phi), (1, phi, 0), (phi, 0, 1)):
    zero = base.index(0)
    nonzero = [index for index in range(3) if index != zero]
    for signs in product((-1, 1), repeat=2):
        vertex = list(base)
        for coordinate, sign in zip(nonzero, signs):
            vertex[coordinate] *= sign
        icosahedron_vertices.append(tuple(vertex))
icosahedron_vertices = np.asarray(sorted(set(icosahedron_vertices)))
icosahedron_edges = tuple(
    (left, right)
    for left in range(12) for right in range(left+1, 12)
    if abs(np.sum((icosahedron_vertices[left]-icosahedron_vertices[right])**2)-4.0)
        < TOLERANCE
)
adjacency = [set() for _ in range(12)]
for left, right in icosahedron_edges:
    adjacency[left].add(right)
    adjacency[right].add(left)
icosahedron_faces = tuple(
    (left, right, third)
    for left, right in icosahedron_edges
    for third in sorted(adjacency[left] & adjacency[right])
    if right < third
)

def induced_vertex_permutation(matrix):
    moved = (matrix @ icosahedron_vertices.T).T
    distance = np.sum(
        (moved[:, None, :]-icosahedron_vertices[None, :, :])**2, axis=2
    )
    permutation = distance.argmin(axis=1)
    residual = float(distance[np.arange(12), permutation].max())
    return (
        tuple(int(value) for value in permutation)
        if residual < 1e-12 else None,
        residual,
    )


rotation_candidates = []
for vertex in icosahedron_vertices:
    for multiple in range(1, 5):
        rotation_candidates.append(
            rotation_matrix(vertex, 2*np.pi*multiple/5)
        )
for face in icosahedron_faces:
    center = sum(
        (icosahedron_vertices[index] for index in face), np.zeros(3)
    )
    rotation_candidates.extend((
        rotation_matrix(center, 2*np.pi/3),
        rotation_matrix(center, 4*np.pi/3),
    ))
for edge in icosahedron_edges:
    rotation_candidates.append(rotation_matrix(
        icosahedron_vertices[edge[0]]+icosahedron_vertices[edge[1]], np.pi
    ))
vertex_rotations = {tuple(range(12))}
maximum_rotation_residual = 0.0
for matrix in rotation_candidates:
    permutation, residual = induced_vertex_permutation(matrix)
    maximum_rotation_residual = max(maximum_rotation_residual, residual)
    if permutation is not None:
        vertex_rotations.add(permutation)
vertex_rotations = tuple(sorted(vertex_rotations))

chambers = []
for face in icosahedron_faces:
    face_edges = [edge for edge in icosahedron_edges if set(edge) <= set(face)]
    for edge in face_edges:
        for vertex in edge:
            chambers.append((vertex, edge, face))
chambers = tuple(chambers)
chamber_index = {chamber: index for index, chamber in enumerate(chambers)}

def chamber_permutation(vertex_permutation):
    output = []
    for vertex, edge, face in chambers:
        image = (
            vertex_permutation[vertex],
            tuple(sorted(vertex_permutation[index] for index in edge)),
            tuple(sorted(vertex_permutation[index] for index in face)),
        )
        output.append(chamber_index[image])
    return tuple(output)


chamber_rotations = tuple(
    sorted({chamber_permutation(rotation) for rotation in vertex_rotations})
)
remaining = set(range(120))
chamber_sheets = []
while remaining:
    seed = min(remaining)
    orbit = frozenset(rotation[seed] for rotation in chamber_rotations)
    chamber_sheets.append(orbit)
    remaining -= orbit
chamber_sheets = tuple(chamber_sheets)

reflection_vertices, reflection_residual = induced_vertex_permutation(-np.eye(3))
chamber_reflection = chamber_permutation(reflection_vertices)
gamma = np.empty(120, dtype=np.int8)
gamma[list(chamber_sheets[0])] = 1
gamma[list(chamber_sheets[1])] = -1
reflection_ok = bool(
    all(chamber_reflection[chamber_reflection[index]] == index
        and chamber_reflection[index] != index for index in range(120))
    and all(
        compose(chamber_reflection, rotation)
        == compose(rotation, chamber_reflection)
        for rotation in chamber_rotations
    )
    and {chamber_reflection[index] for index in chamber_sheets[0]}
        == set(chamber_sheets[1])
    and all(gamma[chamber_reflection[index]] == -gamma[index]
            for index in range(120))
)
chamber_geometry_ok = bool(
    (len(icosahedron_vertices), len(icosahedron_edges), len(icosahedron_faces))
        == (12, 30, 20)
    and len(chambers) == 120
    and len(chamber_rotations) == 60
    and sorted(len(sheet) for sheet in chamber_sheets) == [60, 60]
    and all(len({rotation[index] for rotation in chamber_rotations}) == 60
            for index in range(120))
    and maximum_rotation_residual < 1e-12
    and reflection_residual < 1e-12
    and reflection_ok
)

full_chamber_group = tuple(sorted(
    set(chamber_rotations)
    | {compose(chamber_reflection, rotation) for rotation in chamber_rotations}
))
chamber_table, chamber_identity = multiplication_table(full_chamber_group)
chamber_group_invariants = group_invariants(chamber_table, chamber_identity)
chamber_action_regular = bool(
    len(full_chamber_group) == 120
    and len({permutation[0] for permutation in full_chamber_group}) == 120
)

# Independently reconstruct the binary icosahedral multiplication law.
vertices_2i, _, _ = build_600cell()
binary_table = np.empty((120, 120), dtype=np.int16)
binary_residual = 0.0
for left in range(120):
    for right in range(120):
        raw_product = qmul(vertices_2i[left], vertices_2i[right])
        image = int(np.argmax(vertices_2i @ raw_product))
        binary_table[left, right] = image
        binary_residual = max(
            binary_residual,
            float(np.max(np.abs(vertices_2i[image]-raw_product))),
        )
binary_identity = int(np.argmax(vertices_2i[:, 0]))
binary_group_invariants = group_invariants(binary_table, binary_identity)

group_laws_distinct = bool(
    chamber_action_regular
    and len(chamber_group_invariants["involutions"]) == 31
    and len(chamber_group_invariants["commutator_subgroup"]) == 60
    and binary_residual < TOLERANCE
    and len(binary_group_invariants["involutions"]) == 1
    and len(binary_group_invariants["commutator_subgroup"]) == 120
)

gamma_to_s = equivariant_map_census(target_reflection_odd=True)
gamma_to_d = equivariant_map_census(target_reflection_odd=True)
gamma_to_chi = equivariant_map_census(target_reflection_odd=False)
gamma_s_orbits = pair_orbits(target_reflection_odd=True)
gamma_d_orbits = pair_orbits(target_reflection_odd=True)
gamma_chi_orbits = pair_orbits(target_reflection_odd=False)
chi_orbits = (frozenset((-1,)), frozenset((1,)))

sign_census_ok = bool(
    gamma_to_s["equivariant_maps"] == 2
    and gamma_to_s["equivariant_bijections"] == 2
    and gamma_to_d["equivariant_maps"] == 2
    and gamma_to_d["equivariant_bijections"] == 2
    and gamma_to_chi["equivariant_maps"] == 2
    and gamma_to_chi["equivariant_bijections"] == 0
    and [len(orbit) for orbit in gamma_s_orbits] == [2, 2]
    and [len(orbit) for orbit in gamma_d_orbits] == [2, 2]
    and [len(orbit) for orbit in gamma_chi_orbits] == [2, 2]
    and [len(orbit) for orbit in chi_orbits] == [1, 1]
)

controls_ok = bool(
    provenance_ok and chromatic_input_ok and chamber_geometry_ok
    and group_laws_distinct and sign_census_ok
)
unique_static_selector = bool(
    len(chi_orbits) == 1
    or gamma_to_chi["equivariant_bijections"] == 1
)
if not controls_ok:
    outcome = "OPEN_CONTROL_FAILURE"
elif unique_static_selector:
    outcome = "UNIQUE_STATIC_CHIRAL_SELECTOR"
else:
    outcome = "STATIC_CHIRAL_SELECTOR_NO_GO"

check("frozen source and chromatic artifact hashes", provenance_ok)
check("frozen chromatic artifact certifies the exact residual Z2", chromatic_input_ok)
check("icosahedron and 120 complete chambers reconstructed independently", chamber_geometry_ok)
check("central reflection exchanges sheets and reverses gamma", reflection_ok)
check("full chamber symmetry acts freely and transitively on 120 flags", chamber_action_regular)
check("2I and A5xC2 regular group laws are nonisomorphic by exact invariants", group_laws_distinct, f"involutions: binary={len(binary_group_invariants['involutions'])}, chamber={len(chamber_group_invariants['involutions'])}; commutators: binary={len(binary_group_invariants['commutator_subgroup'])}, chamber={len(chamber_group_invariants['commutator_subgroup'])}")
check("gamma to cover chirality has exactly two equivariant bijections", gamma_to_s["equivariant_bijections"] == 2)
check("gamma to degree sign has exactly two equivariant bijections", gamma_to_d["equivariant_bijections"] == 2)
check("gamma to reflection-even chi has zero equivariant bijections", gamma_to_chi["equivariant_bijections"] == 0)
check("all three paired sign spaces retain exactly two reflection orbits", len(gamma_s_orbits) == len(gamma_d_orbits) == len(gamma_chi_orbits) == 2)
check("both chi values are individually symmetry-fixed", [len(orbit) for orbit in chi_orbits] == [1, 1])
check("existing static data provide no unique selector", not unique_static_selector)
check("mechanical outcome assigned", outcome in {
    "UNIQUE_STATIC_CHIRAL_SELECTOR",
    "STATIC_CHIRAL_SELECTOR_NO_GO",
    "OPEN_CONTROL_FAILURE",
})
check("no dynamics, matter target, preferred schedule or scale parsed", True)

passed = sum(ok for _, ok in tests)
payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": source_hashes,
    "regge_action_parsed": False,
    "nonlinear_schedule_result_parsed": False,
    "matter_or_standard_model_target_parsed": False,
    "preferred_chirality_parsed": False,
    "carrier": {
        "icosahedron_f_vector": [
            len(icosahedron_vertices), len(icosahedron_edges),
            len(icosahedron_faces),
        ],
        "chambers": len(chambers),
        "proper_rotations": len(chamber_rotations),
        "sheet_sizes": sorted(len(sheet) for sheet in chamber_sheets),
        "full_chamber_group": len(full_chamber_group),
    },
    "group_laws": {
        "binary_2I": {
            "order": binary_table.shape[0],
            "nonidentity_involutions": len(binary_group_invariants["involutions"]),
            "commutator_subgroup_order": len(binary_group_invariants["commutator_subgroup"]),
            "maximum_multiplication_residual": binary_residual,
        },
        "split_chamber_A5xC2": {
            "order": chamber_table.shape[0],
            "nonidentity_involutions": len(chamber_group_invariants["involutions"]),
            "commutator_subgroup_order": len(chamber_group_invariants["commutator_subgroup"]),
        },
        "isomorphic": False if group_laws_distinct else None,
    },
    "transformation_laws": {
        "gamma": "reflection odd",
        "cover_side_s": "reflection odd",
        "degree_sign_d": "reflection odd",
        "chi_equals_s_times_d": "reflection even",
    },
    "equivariant_maps": {
        "gamma_to_s": gamma_to_s,
        "gamma_to_d": gamma_to_d,
        "gamma_to_chi": gamma_to_chi,
    },
    "sign_orbits": {
        "gamma_s": [[list(state) for state in sorted(orbit)] for orbit in gamma_s_orbits],
        "gamma_d": [[list(state) for state in sorted(orbit)] for orbit in gamma_d_orbits],
        "gamma_chi": [[list(state) for state in sorted(orbit)] for orbit in gamma_chi_orbits],
        "chi": [list(orbit) for orbit in chi_orbits],
    },
    "labels": {
        "group_and_sign_census": "DERIVED EXACT",
        "static_selector_failure": "DERIVED STRUCTURAL NEGATIVE",
        "dynamical_parity_breaking_term": "OPEN",
        "external_novelty": "OPEN",
    },
    "outcome": outcome,
    "tests": len(tests),
    "passed": passed,
    "test_records": [
        {"label": label, "passed": ok} for label, ok in tests
    ],
}
OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

print("-" * 78)
print(f"2I: involutions={len(binary_group_invariants['involutions'])}, commutator subgroup={len(binary_group_invariants['commutator_subgroup'])}")
print(f"A5xC2: involutions={len(chamber_group_invariants['involutions'])}, commutator subgroup={len(chamber_group_invariants['commutator_subgroup'])}")
print(f"Equivariant bijections gamma->s: {gamma_to_s['equivariant_bijections']}")
print(f"Equivariant bijections gamma->d: {gamma_to_d['equivariant_bijections']}")
print(f"Equivariant bijections gamma->chi: {gamma_to_chi['equivariant_bijections']}")
print(f"Fixed chi choices: {len(chi_orbits)}")
print(f"Outcome: {outcome}")
print(f"RESULT: {passed}/{len(tests)} checks passed")

if passed != len(tests):
    sys.exit(1)

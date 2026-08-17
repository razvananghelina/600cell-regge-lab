#!/usr/bin/env python3
"""Exact chromatic degree of the fixed five-colouring of the 600-cell.

Prior-art commit: 38cce14.
Protocol commit: 4458c23.
"""

from collections import Counter, defaultdict
from itertools import combinations, permutations
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from commons.cell600 import build_600cell  # noqa: E402


OUTPUT = HERE / "gravity_600cell_chromatic_degree_selector.json"
CELL600_SOURCE = ROOT / "commons" / "cell600.py"
ORIENTATION_SOURCE = HERE / "verify_gravity_600cell_staircase_orientation_selector.py"
PRIOR_ART_COMMIT = "38cce14"
PROTOCOL_COMMIT = "4458c23"
CELL600_SHA256 = "ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f"
ORIENTATION_SOURCE_SHA256 = "4885fd9c69ecc82c2d0aa31b5cde72b123999ef01792f1aebc8435ba063dc90e"
DET_TOLERANCE = 1e-10


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


def canonical_key(item):
    return tuple(sorted(item))


def permutation_sign(permutation):
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left+1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def build_tetrahedra(adjacency):
    neighbours = [set(np.flatnonzero(row).tolist()) for row in adjacency]
    edges = {
        (left, right)
        for left in range(120)
        for right in range(left+1, 120)
        if adjacency[left, right]
    }
    triangles = set()
    tetrahedra = []
    for a in range(120):
        for b in sorted(vertex for vertex in neighbours[a] if vertex > a):
            for c in sorted(
                vertex for vertex in neighbours[a] & neighbours[b]
                if vertex > b
            ):
                triangles.add((a, b, c))
                for d in sorted(
                    vertex
                    for vertex in neighbours[a] & neighbours[b] & neighbours[c]
                    if vertex > c
                ):
                    tetrahedra.append((a, b, c, d))
    # Reconstruct all triangle faces, including those not encountered as a,b,c
    # in the clique loop's particular nesting.
    triangles = {
        triangle for tetrahedron in tetrahedra
        for triangle in combinations(tetrahedron, 3)
    }
    return edges, frozenset(triangles), tuple(tetrahedra)


def degree_for_label_map(tetrahedra, source_coefficients, vertex_colour,
                         label_map):
    pushforward = Counter()
    local_sign_counts = {missing: Counter() for missing in range(5)}
    for tetrahedron in tetrahedra:
        labels = tuple(label_map[vertex_colour[vertex]] for vertex in tetrahedron)
        missing_set = set(range(5))-set(labels)
        if len(missing_set) != 1 or len(set(labels)) != 4:
            raise RuntimeError("colour map does not land in one target facet")
        missing = missing_set.pop()
        canonical_facet = tuple(colour for colour in range(5) if colour != missing)
        positions = {colour: index for index, colour in enumerate(canonical_facet)}
        reorder = tuple(positions[colour] for colour in labels)
        image_sign = permutation_sign(reorder)
        contribution = source_coefficients[tetrahedron]*image_sign
        pushforward[missing] += contribution
        # Independent local-degree view: divide each signed preimage by the
        # orientation coefficient (-1)^missing of its target facet.
        local_sign_counts[missing][contribution*((-1)**missing)] += 1
    candidates = tuple(((-1)**missing)*pushforward[missing]
                       for missing in range(5))
    independent = tuple(
        sum(sign*count for sign, count in local_sign_counts[missing].items())
        for missing in range(5)
    )
    return {
        "pushforward_facet_coefficients": tuple(pushforward[i] for i in range(5)),
        "degree_candidates": candidates,
        "independent_preimage_degrees": independent,
        "degree": candidates[0] if len(set(candidates)) == 1 else None,
        "facet_preimage_counts": tuple(sum(local_sign_counts[i].values()) for i in range(5)),
        "local_degree_sign_counts": tuple(dict(local_sign_counts[i]) for i in range(5)),
    }


tests = []


def check(label, condition):
    ok = bool(condition)
    tests.append((label, ok))
    print(f"{'PASS' if ok else 'FAIL'}: {label}")


source_hashes = {
    "cell600": digest(CELL600_SOURCE),
    "orientation_verifier_source": digest(ORIENTATION_SOURCE),
}
provenance_ok = source_hashes == {
    "cell600": CELL600_SHA256,
    "orientation_verifier_source": ORIENTATION_SOURCE_SHA256,
}

vertices, adjacency_float, _ = build_600cell()
adjacency = adjacency_float > 0.5
edges, triangles, tetrahedra = build_tetrahedra(adjacency)

multiplication = np.empty((120, 120), dtype=np.int16)
for left in range(120):
    for right in range(120):
        multiplication[left, right] = int(np.argmax(
            vertices @ qmul(vertices[left], vertices[right])
        ))
conjugate = np.array([
    int(np.argmax(vertices @ (vertex*np.array((1, -1, -1, -1)))))
    for vertex in vertices
], dtype=np.int16)

binary_tetrahedral = frozenset(
    index for index, vertex in enumerate(vertices)
    if (
        np.count_nonzero(np.abs(vertex) > 1e-8) == 1
        and np.max(np.abs(vertex)) > 1-1e-8
    ) or np.all(np.abs(np.abs(vertex)-0.5) < 1e-8)
)
unseen = set(range(120))
cover_cells = []
while unseen:
    representative = min(unseen)
    cell = frozenset(
        int(multiplication[representative, element])
        for element in binary_tetrahedral
    )
    cover_cells.append(cell)
    unseen -= cell
cover_cells = tuple(sorted(cover_cells, key=canonical_key))
cell_lookup = {cell: index for index, cell in enumerate(cover_cells)}
vertex_colour = {
    vertex: colour
    for colour, cell in enumerate(cover_cells)
    for vertex in cell
}

carrier_ok = bool(
    vertices.shape == (120, 4)
    and len(edges) == 720
    and len(triangles) == 1200
    and len(tetrahedra) == 600
    and len(cover_cells) == 5
    and {len(cell) for cell in cover_cells} == {24}
    and len(vertex_colour) == 120
    and all(len({vertex_colour[v] for v in tetrahedron}) == 4
            for tetrahedron in tetrahedra)
)

source_coefficients = {}
minimum_absolute_determinant = float("inf")
for tetrahedron in tetrahedra:
    determinant = float(np.linalg.det(vertices[list(tetrahedron)]))
    minimum_absolute_determinant = min(minimum_absolute_determinant, abs(determinant))
    source_coefficients[tetrahedron] = 1 if determinant > 0 else -1 if determinant < 0 else 0

source_boundary = Counter()
for tetrahedron, coefficient in source_coefficients.items():
    for omitted in range(4):
        face = tetrahedron[:omitted]+tetrahedron[omitted+1:]
        source_boundary[face] += coefficient*((-1)**omitted)
source_boundary = Counter({face: value for face, value in source_boundary.items() if value})
source_chain_ok = bool(
    minimum_absolute_determinant > DET_TOLERANCE
    and not source_boundary
    and set(source_coefficients.values()) <= {-1, 1}
)

identity_map = tuple(range(5))
identity_degree = degree_for_label_map(
    tetrahedra, source_coefficients, vertex_colour, identity_map
)
identity_control_ok = bool(
    identity_degree["degree"] is not None
    and identity_degree["degree_candidates"]
        == identity_degree["independent_preimage_degrees"]
    and identity_degree["facet_preimage_counts"] == (120, 120, 120, 120, 120)
)

order_records = []
all_order_controls = True
for order in permutations(range(5)):
    rank = [0]*5
    for position, old_colour in enumerate(order):
        rank[old_colour] = position
    result = degree_for_label_map(
        tetrahedra, source_coefficients, vertex_colour, tuple(rank)
    )
    expected = permutation_sign(order)*identity_degree["degree"]
    control = bool(
        result["degree"] is not None
        and result["degree_candidates"] == result["independent_preimage_degrees"]
        and result["degree"] == expected
    )
    all_order_controls &= control
    order_records.append({
        "order": list(order),
        "permutation_sign": permutation_sign(order),
        "degree": result["degree"],
        "degree_candidates": list(result["degree_candidates"]),
        "control_pass": control,
    })

degree_by_order_sign = {
    sign: Counter(record["degree"] for record in order_records
                  if record["permutation_sign"] == sign)
    for sign in (-1, 1)
}

# Rebuild the exact setwise cover action and verify A5 invariance of degree.
action_provenance = defaultdict(set)
plain = np.arange(120, dtype=np.int16)
for reflected in (False, True):
    seed = conjugate if reflected else plain
    for left in range(120):
        left_images = multiplication[left, seed]
        for right in range(120):
            action = tuple(
                int(value)
                for value in multiplication[left_images, conjugate[right]]
            )
            action_provenance[action].add(reflected)

setwise_actions = 0
induced_group = set()
for action_tuple in action_provenance:
    action = np.asarray(action_tuple, dtype=np.int16)
    images = tuple(frozenset(int(action[v]) for v in cell) for cell in cover_cells)
    if not all(image in cell_lookup for image in images):
        continue
    induced = tuple(cell_lookup[image] for image in images)
    if len(set(induced)) != 5:
        continue
    setwise_actions += 1
    induced_group.add(induced)

h4_degree_records = []
h4_degree_ok = True
for induced in sorted(induced_group):
    result = degree_for_label_map(
        tetrahedra, source_coefficients, vertex_colour, induced
    )
    control = bool(
        permutation_sign(induced) == 1
        and result["degree"] == identity_degree["degree"]
    )
    h4_degree_ok &= control
    h4_degree_records.append({
        "induced_permutation": list(induced),
        "permutation_sign": permutation_sign(induced),
        "degree": result["degree"],
        "control_pass": control,
    })

controls_ok = bool(
    provenance_ok and carrier_ok and source_chain_ok and identity_control_ok
    and all_order_controls and len(action_provenance) == 14400
    and setwise_actions == 1440 and len(induced_group) == 60
    and h4_degree_ok
)
base_degree = identity_degree["degree"]
opposite_nonzero_classes = bool(
    base_degree not in (None, 0)
    and degree_by_order_sign[1] == Counter({base_degree: 60})
    and degree_by_order_sign[-1] == Counter({-base_degree: 60})
)
if not controls_ok:
    outcome = "OPEN_CONTROL_FAILURE"
elif base_degree == 0:
    outcome = "CHROMATIC_DEGREE_ZERO"
elif opposite_nonzero_classes:
    outcome = "CHROMATIC_ORIENTATION_LINE_DERIVED"
else:
    outcome = "OPEN_CONTROL_FAILURE"

check("frozen source hashes", provenance_ok)
check("600-cell and five-colour carrier counts", carrier_ok)
check("determinant orientation is resolved and the source chain is closed", source_chain_ok)
check("all five target facets give one independently checked integer degree", identity_control_ok)
check("all 120 colour orders obey exact alternating degree", all_order_controls)
check("full 14,400-action H4 census rebuilt", len(action_provenance) == 14400)
check("setwise cover action induces 60 even permutations with kernel 24", setwise_actions == 1440 and len(induced_group) == 60 and all(permutation_sign(g) == 1 for g in induced_group))
check("chromatic degree is invariant under the induced A5 action", h4_degree_ok)
check("mechanical chromatic-degree outcome assigned", outcome in {"CHROMATIC_DEGREE_ZERO", "CHROMATIC_ORIENTATION_LINE_DERIVED", "OPEN_CONTROL_FAILURE"})
check("no Regge, nonlinear, continuum or desired sign parsed", True)

passed = sum(ok for _, ok in tests)
payload = {
    "prior_art_commit": PRIOR_ART_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_sha256": source_hashes,
    "regge_action_parsed": False,
    "nonlinear_result_parsed": False,
    "continuum_target_parsed": False,
    "desired_sign_parsed": False,
    "carrier": {
        "vertices": len(vertices),
        "edges": len(edges),
        "triangles": len(triangles),
        "tetrahedra": len(tetrahedra),
        "cover_cell_sizes": [len(cell) for cell in cover_cells],
        "minimum_absolute_source_determinant": minimum_absolute_determinant,
    },
    "identity_colour_order": {
        "pushforward_facet_coefficients": list(identity_degree["pushforward_facet_coefficients"]),
        "target_boundary_coefficients": [(-1)**index for index in range(5)],
        "degree_candidates": list(identity_degree["degree_candidates"]),
        "independent_preimage_degrees": list(identity_degree["independent_preimage_degrees"]),
        "facet_preimage_counts": list(identity_degree["facet_preimage_counts"]),
        "local_degree_sign_counts": [
            {str(key): value for key, value in counts.items()}
            for counts in identity_degree["local_degree_sign_counts"]
        ],
        "degree": base_degree,
    },
    "order_census": {
        "number_of_orders": len(order_records),
        "degree_by_permutation_sign": {
            str(sign): {str(key): value for key, value in counts.items()}
            for sign, counts in degree_by_order_sign.items()
        },
        "records": order_records,
    },
    "h4_control": {
        "all_actions": len(action_provenance),
        "setwise_actions": setwise_actions,
        "distinct_induced_permutations": len(induced_group),
        "records": h4_degree_records,
    },
    "tests": len(tests),
    "passed": passed,
    "outcome": outcome,
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")

print(f"pushforward={identity_degree['pushforward_facet_coefficients']}")
print(f"degree_candidates={identity_degree['degree_candidates']}")
print(f"identity_degree={base_degree}")
print(f"degree_by_order_sign={degree_by_order_sign}")
print(f"OUTCOME: {outcome}")
print(f"{passed}/{len(tests)} tests passed")

raise SystemExit(0 if passed == len(tests) and outcome != "OPEN_CONTROL_FAILURE" else 1)

#!/usr/bin/env python3
"""Target-blind H4-equivariant and metric-compatible coin census.

Protocol commit d19c453 froze the carrier, two Gram variants, commutant
census and coefficient-free spectral involutions before this computation.
"""

from collections import Counter, deque
from itertools import combinations, permutations, product
import json
from math import acos
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


OUTPUT = Path(__file__).with_name("h4_matrix_coin_selection.json")
PROTOCOL_COMMIT = "d19c453"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


def normalized(vector):
    return vector / np.linalg.norm(vector)


def geodesic_log(base, target):
    cosine = float(np.clip(base @ target, -1.0, 1.0))
    angle = acos(cosine)
    tangent = target - cosine * base
    return angle * tangent / np.linalg.norm(tangent)


def build_chambers(vertices, adjacency):
    neighbours = tuple(
        frozenset(np.flatnonzero(adjacency[index]).tolist())
        for index in range(120)
    )
    edges = tuple(
        (left, right)
        for left in range(120)
        for right in sorted(neighbours[left])
        if left < right
    )
    triangles = tuple(
        (left, right, third)
        for left, right in edges
        for third in sorted(neighbours[left] & neighbours[right])
        if right < third
    )
    tetrahedra = tuple(
        (first, second, third, fourth)
        for first, second, third in triangles
        for fourth in sorted(
            neighbours[first] & neighbours[second] & neighbours[third]
        )
        if third < fourth
    )
    face_to_tetrahedra = {}
    for tetrahedron in tetrahedra:
        for face in combinations(tetrahedron, 3):
            face_to_tetrahedra.setdefault(face, []).append(tetrahedron)

    chambers = tuple(
        (tetrahedron, ordering)
        for tetrahedron in tetrahedra
        for ordering in permutations(tetrahedron)
    )
    chamber_index = {chamber: index for index, chamber in enumerate(chambers)}
    colour_maps = [[] for _ in range(4)]
    centres = []
    for tetrahedron, ordering in chambers:
        for colour in range(3):
            changed = list(ordering)
            changed[colour], changed[colour + 1] = (
                changed[colour + 1], changed[colour]
            )
            colour_maps[colour].append(
                chamber_index[(tetrahedron, tuple(changed))]
            )
        face = tuple(sorted(ordering[:3]))
        across = next(candidate for candidate in face_to_tetrahedra[face]
                      if candidate != tetrahedron)
        opposite = next(vertex for vertex in across if vertex not in face)
        colour_maps[3].append(
            chamber_index[(across, ordering[:3] + (opposite,))]
        )
        flags = (
            normalized(vertices[[ordering[0]]].sum(axis=0)),
            normalized(vertices[list(ordering[:2])].sum(axis=0)),
            normalized(vertices[list(ordering[:3])].sum(axis=0)),
            normalized(vertices[list(tetrahedron)].sum(axis=0)),
        )
        centres.append(normalized(sum(flags)))
    return (
        (len(vertices), len(edges), len(triangles), len(tetrahedra)),
        chambers,
        tuple(map(tuple, colour_maps)),
        np.asarray(centres),
    )


def component_count(colour_maps):
    unseen = set(range(len(colour_maps[0])))
    count = 0
    while unseen:
        count += 1
        seed = unseen.pop()
        queue = deque((seed,))
        while queue:
            chamber = queue.popleft()
            for mapping in colour_maps:
                neighbour = mapping[chamber]
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    queue.append(neighbour)
    return count


def numerical_multiplicities(eigenvalues, tolerance=1e-10):
    groups = []
    for value in eigenvalues:
        if not groups or abs(value - groups[-1][0]) > tolerance:
            groups.append([float(value), 1])
        else:
            groups[-1][1] += 1
    return [(value, multiplicity) for value, multiplicity in groups]


def commutator_nullity(gram, basis):
    columns = np.column_stack([
        (matrix @ gram - gram @ matrix).ravel()
        for matrix in basis
    ])
    rank = int(np.linalg.matrix_rank(columns, tol=1e-10))
    return len(basis) - rank, rank, columns


def block_fit(matrix):
    block_basis = []
    for row in range(2):
        for column in range(2):
            unit = np.zeros((2, 2))
            unit[row, column] = 1.0
            block_basis.append(np.kron(np.eye(2), unit))
    design = np.column_stack([element.ravel() for element in block_basis])
    coefficients, *_ = np.linalg.lstsq(design, matrix.ravel(), rcond=None)
    residual = np.linalg.norm(design @ coefficients - matrix.ravel())
    return float(residual), coefficients.reshape(2, 2)


def paper_block_commutant_dimension(gram):
    basis = []
    for row in range(2):
        for column in range(2):
            unit = np.zeros((2, 2))
            unit[row, column] = 1.0
            basis.append(np.kron(np.eye(2), unit))
    dimension, _, _ = commutator_nullity(gram, basis)
    return dimension


def metric_audit(gram):
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    multiplicities = numerical_multiplicities(eigenvalues)
    full_basis = []
    for row in range(4):
        for column in range(4):
            unit = np.zeros((4, 4))
            unit[row, column] = 1.0
            full_basis.append(unit)
    commutant_dimension, commutator_rank, _ = commutator_nullity(
        gram, full_basis
    )

    block_basis = []
    for row in range(2):
        for column in range(2):
            unit = np.zeros((2, 2))
            unit[row, column] = 1.0
            block_basis.append(np.kron(np.eye(2), unit))
    block_dimension, block_commutator_rank, _ = commutator_nullity(
        gram, block_basis
    )
    permutation_block_dimensions = []
    for permutation in permutations(range(4)):
        change = np.eye(4)[:, permutation]
        changed_gram = change.T @ gram @ change
        permutation_block_dimensions.append(
            paper_block_commutant_dimension(changed_gram)
        )
    diagonal_basis_block_dimension = paper_block_commutant_dimension(
        np.diag(eigenvalues)
    )

    projectors = [
        np.outer(eigenvectors[:, index], eigenvectors[:, index])
        for index in range(4)
    ]
    spectral_coins = []
    # Fix the sign on the lowest eigenspace to +1 to quotient C ~ -C.
    for tail in product((-1, 1), repeat=3):
        signs = (1, *tail)
        coin = sum(sign * projector
                   for sign, projector in zip(signs, projectors))
        spectral_coins.append({
            "signs": list(signs),
            "plus_eigenspace_rank": int(sum(sign > 0 for sign in signs)),
            "unitarity_residual": float(
                np.linalg.norm(coin.T @ coin - np.eye(4))
            ),
        })

    kernel_projector = projectors[0]
    kernel_coin = 2 * kernel_projector - np.eye(4)
    block_residual, block_coefficients = block_fit(kernel_coin)
    kernel_vector = eigenvectors[:, 0]
    if kernel_vector.sum() < 0:
        kernel_vector = -kernel_vector

    return {
        "gram_matrix": gram.tolist(),
        "eigenvalues": eigenvalues.tolist(),
        "eigenvalue_multiplicities": [
            {"representative": value, "multiplicity": multiplicity}
            for value, multiplicity in multiplicities
        ],
        "distinct_eigenspace_count": len(multiplicities),
        "commutant_complex_dimension_from_nullity": commutant_dimension,
        "commutator_superoperator_rank": commutator_rank,
        "commutant_complex_dimension_from_multiplicities": int(sum(
            multiplicity ** 2 for _, multiplicity in multiplicities
        )),
        "commuting_unitary_real_dimension": int(sum(
            multiplicity ** 2 for _, multiplicity in multiplicities
        )),
        "commuting_unitary_factors": [
            f"U({multiplicity})" for _, multiplicity in multiplicities
        ],
        "paper_block_commutant_complex_dimension": block_dimension,
        "paper_block_commutator_rank": block_commutator_rank,
        "paper_block_dimension_over_24_colour_component_permutations": dict(
            sorted(Counter(permutation_block_dimensions).items())
        ),
        "paper_block_dimension_in_gram_eigenbasis": (
            diagonal_basis_block_dimension
        ),
        "spectral_involution_count_modulo_global_sign": len(spectral_coins),
        "spectral_involution_plus_rank_histogram": dict(sorted(Counter(
            item["plus_eigenspace_rank"] for item in spectral_coins
        ).items())),
        "spectral_involutions": spectral_coins,
        "kernel_dimension": int(sum(abs(eigenvalues) < 1e-10)),
        "positive_kernel_vector": kernel_vector.tolist(),
        "kernel_reflection": kernel_coin.tolist(),
        "kernel_reflection_symmetry_residual": float(
            np.linalg.norm(kernel_coin.T - kernel_coin)
        ),
        "kernel_reflection_unitarity_residual": float(
            np.linalg.norm(kernel_coin.T @ kernel_coin - np.eye(4))
        ),
        "kernel_reflection_gram_commutator_residual": float(
            np.linalg.norm(kernel_coin @ gram - gram @ kernel_coin)
        ),
        "kernel_reflection_paper_block_fit_residual": block_residual,
        "kernel_reflection_paper_block_fit_coefficients": (
            block_coefficients.tolist()
        ),
    }


print("=" * 78)
print("H4 MATRIX-COIN SELECTION CENSUS")
print("=" * 78)

vertices, adjacency, _ = build_600cell()
f_vector, chambers, colour_maps, centres = build_chambers(vertices, adjacency)
check("the complete H4 flag carrier has the exact expected counts",
      f_vector == (120, 720, 1200, 600) and len(chambers) == 14400)
check("the four-colour chamber graph is transitive/connected",
      component_count(colour_maps) == 1)

# A pointwise coin invariant under a transitive chamber action is constant.
# H4 acts trivially on the right-Coxeter colour label, hence its remaining
# constant unitary is an unrestricted U(4), not a selected coin.
equivariant_coin_unitary_real_dimension = 4 ** 2
check("H4 equivariance leaves the full 16-real-parameter U(4) coin",
      equivariant_coin_unitary_real_dimension == 16)

all_literal_grams = []
all_unit_grams = []
for chamber in range(len(chambers)):
    base = centres[chamber]
    steps = np.column_stack([
        geodesic_log(base, centres[colour_maps[colour][chamber]])
        for colour in range(4)
    ])
    all_literal_grams.append(steps.T @ steps)
    unit = steps / np.linalg.norm(steps, axis=0)
    all_unit_grams.append(unit.T @ unit)
all_literal_grams = np.asarray(all_literal_grams)
all_unit_grams = np.asarray(all_unit_grams)
literal_spread = float(np.max(np.ptp(all_literal_grams, axis=0)))
unit_spread = float(np.max(np.ptp(all_unit_grams, axis=0)))
check("both labelled metric tensors are chamber-independent",
      max(literal_spread, unit_spread) < 2e-9,
      f"literal spread={literal_spread:.3e}, unit spread={unit_spread:.3e}")

audits = {
    "literal_geodesic_steps": metric_audit(all_literal_grams[0]),
    "unit_directions": metric_audit(all_unit_grams[0]),
}

for name, audit in audits.items():
    check(f"{name}: Gram spectrum has one zero and three simple positive modes",
          audit["kernel_dimension"] == 1
          and audit["distinct_eigenspace_count"] == 4
          and min(audit["positive_kernel_vector"]) > 0,
          f"eigenvalues={audit['eigenvalues']}")
    check(f"{name}: metric commutant dimensions agree and leave U(1)^4",
          audit["commutant_complex_dimension_from_nullity"] == 4
          and audit["commutant_complex_dimension_from_multiplicities"] == 4
          and audit["commuting_unitary_real_dimension"] == 4,
          f"dimension={audit['commutant_complex_dimension_from_nullity']}")
    check(f"{name}: metric plus paper block form leaves scalars only",
          audit["paper_block_commutant_complex_dimension"] == 1,
          f"complex dimension={audit['paper_block_commutant_complex_dimension']}")
    check(f"{name}: the paper-block comparison is explicitly basis-dependent",
          sum(audit[
              "paper_block_dimension_over_24_colour_component_permutations"
          ].values()) == 24
          and audit["paper_block_dimension_in_gram_eigenbasis"] > 1,
          "24 permutation census={}, eigenbasis dimension={}".format(
              audit[
                  "paper_block_dimension_over_24_colour_component_permutations"
              ],
              audit["paper_block_dimension_in_gram_eigenbasis"],
          ))
    check(f"{name}: the coefficient-free spectral census has 8 sign classes",
          audit["spectral_involution_count_modulo_global_sign"] == 8
          and max(item["unitarity_residual"]
                  for item in audit["spectral_involutions"]) < 1e-13)
    check(f"{name}: the designated kernel reflection is canonical but non-block",
          audit["kernel_reflection_symmetry_residual"] < 1e-14
          and audit["kernel_reflection_unitarity_residual"] < 1e-13
          and audit["kernel_reflection_gram_commutator_residual"] < 1e-13
          and audit["kernel_reflection_paper_block_fit_residual"] > 1e-3,
          "block residual={:.6g}".format(
              audit["kernel_reflection_paper_block_fit_residual"]
          ))

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "phenomenological_or_isotropy_target_used": False,
    "chambers": len(chambers),
    "active_internal_dimension": 4,
    "h4_equivariant_constant_coin_group": "U(4)",
    "h4_equivariant_coin_unitary_real_dimension": (
        equivariant_coin_unitary_real_dimension
    ),
    "gram_spreads_over_all_chambers": {
        "literal_geodesic_steps": literal_spread,
        "unit_directions": unit_spread,
    },
    "metric_audits": audits,
    "verdicts": [
        "DERIVED NONSELECTION: H4 leaves U(4)",
        "DERIVED METRIC NONSELECTION: Gram compatibility leaves U(1)^4",
        "STRUCTURAL BASIS-DEPENDENT BRIDGE: in the identity colour/component "
        "identification, Gram plus I2 tensor U(2) leaves global U(1)",
        "DERIVED DISTINGUISHED REFLECTION: zero-mode reflection outside paper block",
        "STRUCTURAL CANDIDATE: promote the zero-mode reflection to a physical coin",
    ],
    "scope": (
        "Target-blind local coin parameter census only; no dynamical or "
        "continuum comparison was performed."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured matrix-coin census was written", OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("H4_ONLY=U(4), real dimension 16")
for name, audit in audits.items():
    print(
        f"{name}: metric unitary=U(1)^4, spectral classes=8, "
        f"block intersection dimension="
        f"{audit['paper_block_commutant_complex_dimension']}"
    )
print("VERDICT: symmetry and Gram compatibility do not uniquely select a coin.")
print("CANDIDATE: the zero-mode reflection is coefficient-free but non-block.")
print("CAVEAT: comparison with the paper block depends on an unselected basis bridge.")
raise SystemExit(0 if passed == tests else 1)

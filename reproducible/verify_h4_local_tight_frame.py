#!/usr/bin/env python3
"""Metric scalar tight-frame feasibility for one barycentric H4 chamber.

Protocol commit 9beb49c froze the two variants, equations, calibration,
numerical standards and scope before this computation was run.
"""

from itertools import combinations, permutations
import json
from math import acos
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
import scipy.linalg as la

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


OUTPUT = Path(__file__).with_name("h4_local_tight_frame.json")
PROTOCOL_COMMIT = "9beb49c"
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

        flag_vertices = (
            normalized(vertices[[ordering[0]]].sum(axis=0)),
            normalized(vertices[list(ordering[:2])].sum(axis=0)),
            normalized(vertices[list(ordering[:3])].sum(axis=0)),
            normalized(vertices[list(tetrahedron)].sum(axis=0)),
        )
        centres.append(normalized(sum(flag_vertices)))
    return (
        (len(vertices), len(edges), len(triangles), len(tetrahedra)),
        chambers,
        tuple(map(tuple, colour_maps)),
        np.asarray(centres),
    )


def double_audit(vectors, projector):
    """Solve drift and simultaneous tight-frame systems in float64."""
    vectors = np.asarray(vectors, dtype=float)
    ambient, count = vectors.shape
    drift_matrix = np.vstack((vectors, np.ones((1, count))))
    drift_rhs = np.r_[np.zeros(ambient), 1.0]
    weights, *_ = np.linalg.lstsq(drift_matrix, drift_rhs, rcond=None)
    drift_singular = np.linalg.svd(drift_matrix, compute_uv=False)
    drift_rank = int(np.linalg.matrix_rank(drift_matrix))
    drift_relative_residual = float(
        np.linalg.norm(drift_matrix @ weights - drift_rhs)
        / np.linalg.norm(drift_rhs)
    )

    moment = sum(
        weights[index] * np.outer(vectors[:, index], vectors[:, index])
        for index in range(count)
    )
    trace = float(np.trace(moment))
    isotropic_scale = trace / 3.0
    isotropy_residual = float(
        np.linalg.norm(moment - isotropic_scale * projector) / trace
    )
    tangent_basis = la.null_space((np.eye(ambient) - projector))
    tangent_moment = tangent_basis.T @ moment @ tangent_basis
    eigenvalues = np.linalg.eigvalsh(tangent_moment)

    rows = ambient + 1 + ambient * ambient
    complete_matrix = np.zeros((rows, count + 1))
    complete_rhs = np.zeros(rows)
    complete_matrix[:ambient, :count] = vectors
    complete_matrix[ambient, :count] = 1.0
    complete_rhs[ambient] = 1.0
    offset = ambient + 1
    for index in range(count):
        complete_matrix[offset:, index] = np.outer(
            vectors[:, index], vectors[:, index]
        ).ravel()
    complete_matrix[offset:, count] = -projector.ravel()
    complete_solution, *_ = np.linalg.lstsq(
        complete_matrix, complete_rhs, rcond=None
    )
    complete_singular = np.linalg.svd(complete_matrix, compute_uv=False)
    complete_residual = float(
        np.linalg.norm(complete_matrix @ complete_solution - complete_rhs)
        / np.linalg.norm(complete_rhs)
    )
    complete_rank = int(np.linalg.matrix_rank(complete_matrix))
    complete_augmented_rank = int(np.linalg.matrix_rank(
        np.column_stack((complete_matrix, complete_rhs))
    ))
    off_diagonal_gram = [
        float(vectors[:, left] @ vectors[:, right])
        for left in range(count)
        for right in range(left + 1, count)
    ]

    return {
        "drift_rank": drift_rank,
        "drift_condition_number": float(
            drift_singular[0] / drift_singular[-1]
        ),
        "drift_relative_residual": drift_relative_residual,
        "drift_weights": weights.tolist(),
        "minimum_drift_weight": float(weights.min()),
        "tangent_moment_eigenvalues": eigenvalues.tolist(),
        "eigenvalue_ratio": float(eigenvalues[-1] / eigenvalues[0]),
        "normalized_traceless_residual": isotropy_residual,
        "complete_rank": complete_rank,
        "complete_augmented_rank": complete_augmented_rank,
        "complete_condition_number": float(
            complete_singular[0] / complete_singular[-1]
        ),
        "complete_relative_residual": complete_residual,
        "complete_solution_weights": complete_solution[:count].tolist(),
        "complete_solution_scale": float(complete_solution[count]),
        "complete_solution_minimum_weight": float(
            complete_solution[:count].min()
        ),
        "off_diagonal_gram_values": off_diagonal_gram,
        "off_diagonal_gram_spread": float(np.ptp(off_diagonal_gram)),
    }


def mp_dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def mp_normalized(vector):
    norm = mp.sqrt(mp_dot(vector, vector))
    return [entry / norm for entry in vector]


def reconstruct_coordinate(value):
    phi = (1 + mp.sqrt(5)) / 2
    positive = [mp.mpf(0), mp.mpf(1), mp.mpf("0.5"),
                phi / 2, 1 / (2 * phi)]
    candidates = positive + [-entry for entry in positive[1:]]
    return min(candidates, key=lambda candidate: abs(float(candidate) - value))


def mp_chamber_centre(chamber, vertices):
    tetrahedron, ordering = chamber
    source = {
        index: [reconstruct_coordinate(float(value)) for value in vertices[index]]
        for index in tetrahedron
    }

    def summed(indices):
        return [sum(source[index][axis] for index in indices) for axis in range(4)]

    flags = [
        mp_normalized(summed(ordering[:1])),
        mp_normalized(summed(ordering[:2])),
        mp_normalized(summed(ordering[:3])),
        mp_normalized(summed(tetrahedron)),
    ]
    return mp_normalized([sum(flag[axis] for flag in flags) for axis in range(4)])


def mp_geodesic_log(base, target):
    cosine = mp_dot(base, target)
    angle = mp.acos(cosine)
    tangent = [target[index] - cosine * base[index] for index in range(4)]
    tangent = mp_normalized(tangent)
    return [angle * entry for entry in tangent]


def mp_audit(vectors, projector):
    """High-precision overdetermined QR cross-check."""
    ambient = len(vectors)
    count = len(vectors[0])
    drift_matrix = mp.matrix(ambient + 1, count)
    drift_rhs = mp.matrix(ambient + 1, 1)
    for row in range(ambient):
        for column in range(count):
            drift_matrix[row, column] = vectors[row][column]
    for column in range(count):
        drift_matrix[ambient, column] = 1
    drift_rhs[ambient] = 1
    weights, _ = mp.qr_solve(drift_matrix, drift_rhs)
    drift_residual = mp.norm(drift_matrix * weights - drift_rhs)

    moment = mp.matrix(ambient, ambient)
    for row in range(ambient):
        for column in range(ambient):
            moment[row, column] = sum(
                weights[index] * vectors[row][index] * vectors[column][index]
                for index in range(count)
            )
    scale = sum(moment[index, index] for index in range(ambient)) / 3
    isotropy_residual = mp.norm(moment - scale * projector) / (
        3 * abs(scale)
    )

    rows = ambient + 1 + ambient * ambient
    complete_matrix = mp.matrix(rows, count + 1)
    complete_rhs = mp.matrix(rows, 1)
    for row in range(ambient):
        for column in range(count):
            complete_matrix[row, column] = vectors[row][column]
    for column in range(count):
        complete_matrix[ambient, column] = 1
    complete_rhs[ambient] = 1
    offset = ambient + 1
    for row in range(ambient):
        for column in range(ambient):
            flat = offset + row * ambient + column
            for index in range(count):
                complete_matrix[flat, index] = (
                    vectors[row][index] * vectors[column][index]
                )
            complete_matrix[flat, count] = -projector[row, column]
    complete_solution, _ = mp.qr_solve(complete_matrix, complete_rhs)
    complete_residual = mp.norm(
        complete_matrix * complete_solution - complete_rhs
    )
    off_diagonal_gram = [
        sum(vectors[row][left] * vectors[row][right]
            for row in range(ambient))
        for left in range(count)
        for right in range(left + 1, count)
    ]
    return {
        "drift_weights": [mp.nstr(value, 50) for value in weights],
        "minimum_drift_weight": mp.nstr(min(weights), 50),
        "drift_relative_residual": mp.nstr(drift_residual, 20),
        "normalized_traceless_residual": mp.nstr(isotropy_residual, 50),
        "complete_relative_residual": mp.nstr(complete_residual, 50),
        "complete_solution_weights": [
            mp.nstr(complete_solution[index], 50) for index in range(count)
        ],
        "complete_solution_scale": mp.nstr(complete_solution[count], 50),
        "off_diagonal_gram_values": [
            mp.nstr(value, 50) for value in off_diagonal_gram
        ],
        "off_diagonal_gram_spread": mp.nstr(
            max(off_diagonal_gram) - min(off_diagonal_gram), 50
        ),
    }


print("=" * 78)
print("H4 LOCAL SCALAR TIGHT-FRAME FEASIBILITY")
print("=" * 78)

vertices, adjacency, _ = build_600cell()
f_vector, chambers, colour_maps, centres = build_chambers(vertices, adjacency)
check("the 600-cell and chamber counts are exact",
      f_vector == (120, 720, 1200, 600) and len(chambers) == 14400)
check("all four colour maps are fixed-point-free involutions",
      all(
          all(mapping[mapping[index]] == index and mapping[index] != index
              for index in range(len(chambers)))
          for mapping in colour_maps
      ))

# Compute all local Gram matrices to ensure chamber 0 is not selected for its
# numerical answer.
all_gram = []
for chamber in range(len(chambers)):
    base = centres[chamber]
    steps = np.column_stack([
        geodesic_log(base, centres[colour_maps[colour][chamber]])
        for colour in range(4)
    ])
    all_gram.append(steps.T @ steps)
all_gram = np.asarray(all_gram)
gram_spread = float(np.max(np.ptp(all_gram, axis=0)))
check("the labelled local step Gram matrix is chamber-independent",
      gram_spread < 2e-9, f"maximum spread={gram_spread:.3e}")

base = centres[0]
literal_steps = np.column_stack([
    geodesic_log(base, centres[colour_maps[colour][0]])
    for colour in range(4)
])
lengths = np.linalg.norm(literal_steps, axis=0)
unit_steps = literal_steps / lengths
projector = np.eye(4) - np.outer(base, base)

# Known-answer calibration in an ordinary three-dimensional tangent space.
control = np.asarray((
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
), dtype=float).T / np.sqrt(3)
control_audit = double_audit(control, np.eye(3))
check("regular-tetrahedron calibration recovers the unique quarter weights",
      control_audit["drift_rank"] == 4
      and max(abs(np.asarray(control_audit["drift_weights"]) - 0.25)) < 2e-15
      and control_audit["complete_relative_residual"] < 2e-15
      and control_audit["normalized_traceless_residual"] < 2e-15)

double_results = {
    "literal_geodesic_steps": double_audit(literal_steps, projector),
    "unit_directions": double_audit(unit_steps, projector),
}

# High-precision reconstruction from the exact coordinate alphabet.
mp.mp.dps = 80
mp_centres = [mp_chamber_centre(chambers[index], vertices) for index in [
    0, *(colour_maps[colour][0] for colour in range(4))
]]
mp_base = mp_centres[0]
mp_literal_columns = [
    mp_geodesic_log(mp_base, mp_centres[colour + 1])
    for colour in range(4)
]
mp_literal = [
    [mp_literal_columns[column][row] for column in range(4)]
    for row in range(4)
]
mp_lengths = [mp.sqrt(sum(entry * entry for entry in column))
              for column in mp_literal_columns]
mp_unit = [
    [mp_literal[row][column] / mp_lengths[column] for column in range(4)]
    for row in range(4)
]
mp_projector = mp.eye(4) - mp.matrix(mp_base) * mp.matrix(mp_base).T
high_results = {
    "literal_geodesic_steps": mp_audit(mp_literal, mp_projector),
    "unit_directions": mp_audit(mp_unit, mp_projector),
}

mp_control_scale = mp.sqrt(3)
mp_control = [
    [mp.mpf(value) / mp_control_scale for value in row]
    for row in (
        (1, 1, -1, -1),
        (1, -1, 1, -1),
        (1, -1, -1, 1),
    )
]
mp_control_audit = mp_audit(mp_control, mp.eye(3))
check("the 80-digit solver independently passes the tetrahedron control",
      float(mp_control_audit["complete_relative_residual"]) < 1e-70
      and max(abs(float(value) - 0.25)
              for value in mp_control_audit["drift_weights"]) < 1e-70)

mp_literal_float = np.asarray(mp_literal, dtype=float)
reconstruction_error = float(np.max(abs(mp_literal_float - literal_steps)))
check("80-digit exact-alphabet reconstruction agrees with float geometry",
      reconstruction_error < 2e-9,
      f"maximum step-coordinate discrepancy={reconstruction_error:.3e}")

for name in ("literal_geodesic_steps", "unit_directions"):
    low = double_results[name]
    high = high_results[name]
    check(f"{name}: drift solution is unique and strictly positive",
          low["drift_rank"] == 4
          and low["minimum_drift_weight"] > 1e-12
          and float(high["minimum_drift_weight"]) > 1e-12,
          f"p={low['drift_weights']}")
    low_complete = low["complete_relative_residual"]
    high_complete = float(high["complete_relative_residual"])
    stable_classification = (
        (low_complete < 1e-10 and high_complete < 1e-30)
        or (low_complete > 1e-12 and high_complete > 1e-12)
    )
    check(f"{name}: simultaneous-system verdict is precision-stable",
          stable_classification,
          f"residual64={low_complete:.6e}, residual80={high_complete:.6e}")

check("both unconstrained simultaneous systems are rank-inconsistent",
      all(
          result["complete_rank"] == 5
          and result["complete_augmented_rank"] == 6
          for result in double_results.values()
      ),
      "; ".join(
          f"{name}: rank(A)={result['complete_rank']}, "
          f"rank([A|b])={result['complete_augmented_rank']}"
          for name, result in double_results.items()
      ))

# With four spanning vectors and positive p, zero drift plus a tight-frame
# second moment implies all six off-diagonal entries v_i.v_j are the same.
# Proof: W=V diag(sqrt(p)) obeys WW^T=cI and ker(W)=span(sqrt(p)), hence
# W^T W=c(I-sqrt(p)sqrt(p)^T).  This is an independent geometric obstruction.
check("both H4 variants violate the four-vector simplex Gram necessity",
      all(result["off_diagonal_gram_spread"] > 1e-3
          for result in double_results.values()),
      "; ".join(
          f"{name}: spread={result['off_diagonal_gram_spread']:.6g}"
          for name, result in double_results.items()
      ))

# A stay-put mass q scales both nonzero moments and cannot change normalized
# tensor shape.  Check the identity numerically for a non-special q.
q = 0.37
reference = double_results["literal_geodesic_steps"]
weights = np.asarray(reference["drift_weights"])
moment = sum(weights[index] * np.outer(literal_steps[:, index],
                                       literal_steps[:, index])
             for index in range(4))
scaled = (1 - q) * moment
shape = np.linalg.norm(moment - np.trace(moment) / 3 * projector) / np.trace(moment)
scaled_shape = np.linalg.norm(
    scaled - np.trace(scaled) / 3 * projector
) / np.trace(scaled)
check("adding a stay-put probability leaves normalized anisotropy unchanged",
      abs(shape - scaled_shape) < 2e-15)


def classify(high):
    residual = float(high["complete_relative_residual"])
    minimum = min(float(value) for value in high["complete_solution_weights"])
    scale = float(high["complete_solution_scale"])
    if residual < 1e-30 and minimum > 1e-12 and scale > 0:
        return "positive_tight_frame"
    if residual > 1e-12:
        return "incompatible"
    return "open_numerical"


classifications = {name: classify(result) for name, result in high_results.items()}
if classifications["literal_geodesic_steps"] == "positive_tight_frame":
    verdict = "DERIVED SCALAR METRIC TIGHT FRAME"
elif classifications["unit_directions"] == "positive_tight_frame":
    verdict = "STRUCTURAL DIRECTION-ONLY TIGHT FRAME"
elif all(value == "incompatible" for value in classifications.values()):
    verdict = "DERIVED NUMERICAL SCALAR NO-GO AT FIRST SCALE"
else:
    verdict = "OPEN NUMERICAL"

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "phenomenological_target_used": False,
    "attempt_count_N": 2,
    "variants": ["literal_geodesic_steps", "unit_directions"],
    "reference_chamber": 0,
    "colour_geodesic_step_lengths": lengths.tolist(),
    "gram_spread_over_all_chambers": gram_spread,
    "exact_alphabet_reconstruction_max_error": reconstruction_error,
    "double_precision": double_results,
    "high_precision_control": mp_control_audit,
    "high_precision_80_decimal": high_results,
    "classifications": classifications,
    "unconstrained_simultaneous_compatible": {
        name: bool(result["complete_rank"]
                   == result["complete_augmented_rank"])
        for name, result in double_results.items()
    },
    "positive_simultaneous_solution": {
        name: classification == "positive_tight_frame"
        for name, classification in classifications.items()
    },
    "dwell_probability_shape_identity_check_q": q,
    "verdict": verdict,
    "scope": (
        "Scalar positive weights on the four fixed local H4 directions at "
        "the first barycentric scale; not a general matrix-coin no-go."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured scalar tight-frame certificate was written", OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
for name in ("literal_geodesic_steps", "unit_directions"):
    audit = double_results[name]
    print(
        f"{name}: p={audit['drift_weights']} "
        f"R={audit['eigenvalue_ratio']:.9f} "
        f"A={audit['normalized_traceless_residual']:.9f} "
        f"complete_residual={audit['complete_relative_residual']:.9e}"
    )
print(f"VERDICT: {verdict}")
print("SCOPE: scalar first-scale gate only; matrix-valued tetrads remain open.")
raise SystemExit(0 if passed == tests else 1)

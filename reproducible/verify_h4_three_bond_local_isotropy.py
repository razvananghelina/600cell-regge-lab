#!/usr/bin/env python3
"""One-period geodesic isotropy gate for the H4 three-bond walk.

Protocol commit 6811eee froze the chamber representatives, fixed coin,
maximally mixed local input, all six schedules and numerical thresholds.
"""

from collections import defaultdict
from itertools import combinations, permutations
import json
from math import acos, pi, sqrt
from pathlib import Path
import sys

import numpy as np
import scipy.linalg as la

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from commons import build_600cell


OUTPUT = Path(__file__).with_name("h4_three_bond_local_isotropy.json")
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
    if angle < 1e-14:
        return np.zeros(4)
    tangent = target - cosine * base
    return angle * tangent / np.linalg.norm(tangent)


def moment_audit(probabilities, centres, initial_chamber):
    base = centres[initial_chamber]
    tangent_basis = la.null_space(base.reshape(1, 4))
    vectors = {
        chamber: geodesic_log(base, centres[chamber])
        for chamber in probabilities
    }
    mean = sum(
        probabilities[chamber] * vector
        for chamber, vector in vectors.items()
    )
    covariance = np.zeros((4, 4))
    raw_second_moment = np.zeros((4, 4))
    for chamber, vector in vectors.items():
        probability = probabilities[chamber]
        raw_second_moment += probability * np.outer(vector, vector)
        centred = vector - mean
        covariance += probability * np.outer(centred, centred)
    tangent_covariance = tangent_basis.T @ covariance @ tangent_basis
    tangent_raw = tangent_basis.T @ raw_second_moment @ tangent_basis
    eigenvalues = np.linalg.eigvalsh(tangent_covariance)
    trace = float(np.trace(tangent_covariance))
    isotropic = np.eye(3) * trace / 3.0
    residual = float(np.linalg.norm(tangent_covariance - isotropic) / trace)
    ratio = float(eigenvalues[-1] / eigenvalues[0])
    return {
        "probability_sum": float(sum(probabilities.values())),
        "support_chambers": len(probabilities),
        "mean_vector": mean.tolist(),
        "mean_norm": float(np.linalg.norm(mean)),
        "tangent_covariance": tangent_covariance.tolist(),
        "tangent_raw_second_moment": tangent_raw.tolist(),
        "tangent_covariance_eigenvalues": eigenvalues.tolist(),
        "eigenvalue_ratio": ratio,
        "normalized_traceless_residual": residual,
        "isotropic_gate": bool(
            np.linalg.norm(mean) < 1e-12
            and abs(ratio - 1.0) < 1e-10
            and residual < 1e-10
        ),
    }


print("=" * 78)
print("ONE-PERIOD LOCAL ISOTROPY OF THE H4 THREE-BOND WALK")
print("=" * 78)

# -------------------------------------------------------------------------
# 600-cell complex, complete flags and their geometric chamber centres.
# -------------------------------------------------------------------------
vertices, adjacency, _ = build_600cell()
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
check("the coarse f-vector is exact",
      (len(vertices), len(edges), len(triangles), len(tetrahedra))
      == (120, 720, 1200, 600))

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

colour_maps = tuple(map(tuple, colour_maps))
centres = np.asarray(centres)
check("14,400 normalized flag-chamber representatives are constructed",
      centres.shape == (14400, 4)
      and np.max(abs(np.linalg.norm(centres, axis=1) - 1.0)) < 2e-15)

# Colour edge lengths must be constant by H4 chamber transitivity.  Check the
# numerical realization rather than assuming the symmetry was implemented.
colour_distances = []
colour_distance_spreads = []
for colour, mapping in enumerate(colour_maps):
    cosines = np.sum(centres * centres[np.asarray(mapping)], axis=1)
    distances = np.arccos(np.clip(cosines, -1.0, 1.0))
    colour_distances.append(float(distances.mean()))
    colour_distance_spreads.append(float(np.ptp(distances)))
check("each rank colour has one constant geodesic step length",
      max(colour_distance_spreads) < 2e-9,
      f"lengths={colour_distances}")

# -------------------------------------------------------------------------
# Fixed published coin and exact sparse-state evolution.
# -------------------------------------------------------------------------
sigma_x = np.array(((0, 1), (1, 0)), dtype=complex)
sigma_z = np.array(((1, 0), (0, -1)), dtype=complex)
identity2 = np.eye(2, dtype=complex)
rotation_x = np.cos(pi / 4) * identity2 - 1j * np.sin(pi / 4) * sigma_x
rotation_z = np.cos(pi / 4) * identity2 - 1j * np.sin(pi / 4) * sigma_z
coin = np.exp(1j * pi / 3) * rotation_z @ rotation_x
coin_hat = np.kron(identity2, coin)
check("the fixed directional coin is unitary and has two dense 2x2 blocks",
      np.max(abs(coin_hat.conj().T @ coin_hat - np.eye(4))) < 2e-15
      and np.count_nonzero(abs(coin_hat) > 1e-14) == 8)


def translate(amplitudes, bond):
    first_colour, second_colour = bond
    result = defaultdict(complex)
    for (chamber, component), amplitude in amplitudes.items():
        if component == 0:
            target = (colour_maps[second_colour][chamber], 2)
        elif component == 1:
            target = (chamber, 3)
        elif component == 2:
            target = (chamber, 0)
        else:
            target = (colour_maps[first_colour][chamber], 1)
        result[target] += amplitude
    return dict(result)


def apply_coin(amplitudes):
    result = defaultdict(complex)
    for (chamber, source_component), amplitude in amplitudes.items():
        for output_component in range(4):
            coefficient = coin_hat[output_component, source_component]
            if abs(coefficient) > 1e-14:
                result[(chamber, output_component)] += coefficient * amplitude
    return {
        state: amplitude for state, amplitude in result.items()
        if abs(amplitude) > 1e-14
    }


def evolve_basis(initial_chamber, initial_component, schedule):
    amplitudes = {(initial_chamber, initial_component): 1.0 + 0.0j}
    for bond in schedule:
        amplitudes = apply_coin(translate(amplitudes, bond))
    return amplitudes


# Known-answer calibration: four regular-tetrahedron directions are a tight
# frame in R3 with covariance I/3.
control_directions = np.asarray((
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
), dtype=float) / sqrt(3)
control_mean = control_directions.mean(axis=0)
control_covariance = sum(
    np.outer(direction, direction) for direction in control_directions
) / 4
control_eigenvalues = np.linalg.eigvalsh(control_covariance)
control_ratio = float(control_eigenvalues[-1] / control_eigenvalues[0])
control_residual = float(
    np.linalg.norm(control_covariance - np.eye(3) / 3)
    / np.trace(control_covariance)
)
check("regular-tetrahedron control returns exact numerical isotropy",
      np.linalg.norm(control_mean) < 1e-15
      and abs(control_ratio - 1.0) < 1e-15
      and control_residual < 1e-15,
      f"R={control_ratio}, A={control_residual}")

bonds = ((0, 1), (1, 2), (2, 3))
schedules = tuple(permutations(bonds))
initial_chamber = 0
schedule_audits = []
for schedule in schedules:
    probabilities = defaultdict(float)
    basis_norms = []
    for initial_component in range(4):
        amplitudes = evolve_basis(
            initial_chamber, initial_component, schedule
        )
        norm = sum(abs(amplitude) ** 2 for amplitude in amplitudes.values())
        basis_norms.append(float(norm))
        for (chamber, _component), amplitude in amplitudes.items():
            probabilities[chamber] += abs(amplitude) ** 2 / 4.0
    audit = moment_audit(probabilities, centres, initial_chamber)
    audit["schedule"] = [list(bond) for bond in schedule]
    audit["basis_state_norms"] = basis_norms
    schedule_audits.append(audit)

check("all 24 evolved basis-state norms and six mixed probabilities equal one",
      all(
          max(abs(value - 1.0) for value in audit["basis_state_norms"])
          < 2e-14
          and abs(audit["probability_sum"] - 1.0) < 2e-14
          for audit in schedule_audits
      ))

designated = bonds
designated_audit = next(
    audit for audit in schedule_audits
    if tuple(map(tuple, audit["schedule"])) == designated
)
isotropy_hits = sum(bool(audit["isotropic_gate"])
                    for audit in schedule_audits)
check("the designated schedule fails the frozen one-period isotropy gate",
      not designated_audit["isotropic_gate"],
      "mean={:.6g}, R={:.6g}, A={:.6g}".format(
          designated_audit["mean_norm"],
          designated_audit["eigenvalue_ratio"],
          designated_audit["normalized_traceless_residual"],
      ))
check("all six schedules fail one-period isotropy without order selection",
      isotropy_hits == 0,
      f"isotropy hits={isotropy_hits}/6")

payload = {
    "protocol_commit": "6811eee",
    "phenomenological_target_used": False,
    "chambers": len(chambers),
    "initial_chamber": initial_chamber,
    "initial_internal_state": "maximally mixed I4/4",
    "coin": "exp(i*pi/3) Rz(pi/2) Rx(pi/2), C_hat=I2 tensor C",
    "colour_geodesic_step_lengths": colour_distances,
    "colour_geodesic_step_length_spreads": colour_distance_spreads,
    "control": {
        "mean_norm": float(np.linalg.norm(control_mean)),
        "eigenvalues": control_eigenvalues.tolist(),
        "eigenvalue_ratio": control_ratio,
        "normalized_traceless_residual": control_residual,
    },
    "schedule_count_N": len(schedules),
    "isotropy_hit_fraction": [isotropy_hits, len(schedules)],
    "designated_schedule": [list(bond) for bond in designated],
    "designated_audit": designated_audit,
    "schedule_audits": schedule_audits,
    "verdict": (
        "ONE-PERIOD ISOTROPY" if designated_audit["isotropic_gate"] else
        "DERIVED NUMERICAL ONE-PERIOD LOCAL ANISOTROPY"
    ),
    "scope": (
        "One period at the first barycentric scale; failure is not a "
        "multistep or refinement no-go."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured local-isotropy certificate was written", OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print(f"ONE_PERIOD_ISOTROPY_HITS={isotropy_hits}/6")
print("DESIGNATED mean={:.9f} R={:.9f} A={:.9f}".format(
    designated_audit["mean_norm"],
    designated_audit["eigenvalue_ratio"],
    designated_audit["normalized_traceless_residual"],
))
print("SCOPE: first barycentric scale and one period only.")
raise SystemExit(0 if passed == tests else 1)

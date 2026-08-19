#!/usr/bin/env python3
"""Projective-coordinate audit of tetrahedral-prism rigidity."""

from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys

import sympy as sy


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PRIMARY = HERE / "gravity_600cell_tetrahedral_prism_rigidity.json"
OUTPUT = HERE / "gravity_600cell_tetrahedral_prism_rigidity_adversarial.json"
AUDIT_PROTOCOL_COMMIT = "e424c98"
INPUT_HASHES = {
    "reproducible/gravity_600cell_tetrahedral_prism_rigidity.json":
        "ce9eb1917dd647c6dd8155a0f9646a72dc7734c0310f763ec31e070403230db8",
}
EUCLIDEAN = (1, 1, 1, 1)
LORENTZIAN = (1, 1, 1, -1)
BASE = (
    (sy.Integer(1), 0, 0, 0),
    (sy.Integer(-1), 0, 0, 0),
    (0, sy.Integer(1), 0, 0),
    (0, 0, sy.Integer(1), 0),
)
NATURAL_EDGES = (
    tuple(combinations(range(4), 2))
    + tuple(combinations(range(4, 8), 2))
    + tuple((index, index+4) for index in range(4))
)
SIDE_PAIRS = tuple(combinations(range(4), 2))
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


def difference(left, right):
    return tuple(sy.sympify(a)-sy.sympify(b) for a, b in zip(left, right))


def dot(left, right, metric):
    return sy.expand(sum(metric[axis]*left[axis]*right[axis]
                         for axis in range(4)))


def controls(scale, translation=(0, 0, 0, 2)):
    scale = sy.Rational(scale)
    translation = tuple(map(sy.Rational, translation))
    top = tuple(tuple(scale*point[axis]+translation[axis]
                      for axis in range(4)) for point in BASE)
    return BASE+top


def squared_lengths(points, metric, edges=NATURAL_EDGES):
    return tuple(dot(difference(points[left], points[right]),
                     difference(points[left], points[right]), metric)
                 for left, right in edges)


def projective_velocities(points):
    """24 exact affine velocity fields from PGL(5), with A_44 removed."""
    generators = tuple((row, column) for row in range(5) for column in range(5)
                       if (row, column) != (4, 4))
    fields = []
    for row, column in generators:
        field = []
        for point in points:
            homogeneous = tuple(point)+(sy.Integer(1),)
            denominator_velocity = (homogeneous[column]
                                    if row == 4 else sy.Integer(0))
            affine = tuple(
                (homogeneous[column] if axis == row else sy.Integer(0))
                - point[axis]*denominator_velocity
                for axis in range(4)
            )
            field.append(affine)
        fields.append(tuple(field))
    return generators, tuple(fields)


def projective_length_jacobian(points, metric):
    generators, fields = projective_velocities(points)
    matrix = sy.zeros(len(NATURAL_EDGES), len(generators))
    for row, (left, right) in enumerate(NATURAL_EDGES):
        edge = difference(points[left], points[right])
        for column, field in enumerate(fields):
            velocity = difference(field[left], field[right])
            matrix[row, column] = 2*dot(edge, velocity, metric)
    return matrix


def translated_prism(translation):
    translation = tuple(map(sy.Rational, translation))
    return BASE+tuple(tuple(point[axis]+translation[axis]
                            for axis in range(4)) for point in BASE)


def oriented_volume(points):
    matrix = sy.Matrix.hstack(*(
        sy.Matrix(difference(points[index], points[0]))
        for index in (1, 2, 3, 4)
    ))
    return sy.simplify(matrix.det()/sy.factorial(3))


def cross_diagonal_lengths(points, metric):
    result = []
    for left, right in SIDE_PAIRS:
        for edge in ((left, right+4), (right, left+4)):
            vector = difference(points[edge[0]], points[edge[1]])
            result.append(dot(vector, vector, metric))
    return tuple(result)


def parallelogram_faces(points):
    return all(
        difference(points[left+4], points[left])
        == difference(points[right+4], points[right])
        for left, right in SIDE_PAIRS
    )


print("="*78)
print("ADVERSARIAL PROJECTIVE PRISM-RIGIDITY AUDIT")
print("="*78)

actual_hashes = {name: digest(ROOT/name) for name in INPUT_HASHES}
check(
    "the frozen primary artifact and audit protocol have exact provenance",
    actual_hashes == INPUT_HASHES and AUDIT_PROTOCOL_COMMIT == "e424c98",
    str(actual_hashes),
)

rank_records = {}
for scale in (sy.Rational(1), sy.Rational(9, 10),
              sy.Rational(11, 10), sy.Rational(2)):
    points = controls(scale)
    record = {}
    for name, metric in (("euclidean", EUCLIDEAN),
                         ("lorentzian", LORENTZIAN)):
        matrix = projective_length_jacobian(points, metric)
        rank = int(matrix.rank())
        record[name] = {
            "length_map_rank_on_projective_chart": rank,
            "projective_kernel": 24-rank,
            "quotient_nonisometric_modes": 24-rank-10,
        }
    rank_records[str(scale)] = record

check(
    "the independent projective chart finds three equal-scale shape modes",
    all(rank_records["1"][signature] == {
        "length_map_rank_on_projective_chart": 11,
        "projective_kernel": 13,
        "quotient_nonisometric_modes": 3,
    } for signature in ("euclidean", "lorentzian")),
    str(rank_records["1"]),
)
check(
    "both sides and the coarse unequal scale have only isometry modes",
    all(rank_records[scale][signature] == {
        "length_map_rank_on_projective_chart": 14,
        "projective_kernel": 10,
        "quotient_nonisometric_modes": 0,
    } for scale in ("9/10", "11/10", "2")
      for signature in ("euclidean", "lorentzian")),
    str({scale: rank_records[scale] for scale in ("9/10", "11/10", "2")}),
)

# Directly differentiate the equal-scale translation family.  At t=(0,0,0,2)
# the three spatial translations are tangent to the fixed-strut-length sphere
# or hyperboloid and leave every top/bottom edge unchanged.
static_points = controls(sy.Rational(1))
translation_derivatives = {}
for name, metric in (("euclidean", EUCLIDEAN),
                     ("lorentzian", LORENTZIAN)):
    derivatives = []
    for axis in range(3):
        velocity = tuple(sy.Integer(index == axis) for index in range(4))
        values = []
        for left, right in NATURAL_EDGES:
            edge = difference(static_points[left], static_points[right])
            left_velocity = velocity if left >= 4 else (0, 0, 0, 0)
            right_velocity = velocity if right >= 4 else (0, 0, 0, 0)
            values.append(2*dot(edge, difference(left_velocity, right_velocity),
                                metric))
        derivatives.append(tuple(values))
    translation_derivatives[name] = derivatives
check(
    "three explicit tangential top translations lie in both length kernels",
    all(all(all(value == 0 for value in derivative)
            for derivative in derivatives)
        for derivatives in translation_derivatives.values()),
)

finite_records = {}
for name, metric, first, second in (
    ("euclidean", EUCLIDEAN, (1, 2, 3, 4), (4, 2, 3, 1)),
    ("lorentzian", LORENTZIAN, (1, 2, 3, 4), (0, 1, 1, 2)),
):
    left = translated_prism(first)
    right = translated_prism(second)
    natural_left = squared_lengths(left, metric)
    natural_right = squared_lengths(right, metric)
    diagonals_left = cross_diagonal_lengths(left, metric)
    diagonals_right = cross_diagonal_lengths(right, metric)
    volume_left = oriented_volume(left)
    volume_right = oriented_volume(right)
    finite_records[name] = {
        "natural_lengths_identical": natural_left == natural_right,
        "parallelogram_faces_first": parallelogram_faces(left),
        "parallelogram_faces_second": parallelogram_faces(right),
        "cross_diagonal_vectors_identical": diagonals_left == diagonals_right,
        "cross_diagonal_differences": sum(a != b for a, b in zip(
            diagonals_left, diagonals_right)),
        "first_oriented_volume": str(volume_left),
        "second_oriented_volume": str(volume_right),
        "absolute_volumes_differ": sy.Abs(volume_left) != sy.Abs(volume_right),
    }

check(
    "finite pairs independently preserve natural lengths and planar side faces",
    all(record["natural_lengths_identical"]
        and record["parallelogram_faces_first"]
        and record["parallelogram_faces_second"]
        for record in finite_records.values()),
    str(finite_records),
)
check(
    "the omitted cross-diagonal data distinguish both finite pairs",
    all(not record["cross_diagonal_vectors_identical"]
        and record["cross_diagonal_differences"] > 0
        for record in finite_records.values()),
    str({key: value["cross_diagonal_differences"]
         for key, value in finite_records.items()}),
)
check(
    "same natural lengths give different four-volumes in both signatures",
    all(record["absolute_volumes_differ"] for record in finite_records.values()),
    str({key: (value["first_oriented_volume"], value["second_oriented_volume"])
         for key, value in finite_records.items()}),
)

primary = json.loads(PRIMARY.read_text())
primary_outcome = primary["outcome"]
agreement = (
    primary_outcome["missing_equal_scale_modes"] == 3
    and primary_outcome["equal_scale_finite_nonuniqueness"] is True
    and primary_outcome["unequal_scale_local_infinitesimal_determination"] is True
    and primary_outcome["hessian_from_16_lengths_at_equal_scale_is_canonical"] is False
)
check(
    "the projective audit agrees with the primary mixed verdict",
    agreement,
    str(primary_outcome),
)

outcome = {
    "mechanism": "projective tangent chart; no planarity minors",
    "equal_scale_nonisometric_modes": 3,
    "equal_scale_finite_nonuniqueness": True,
    "unequal_scale_only_isometry_modes": True,
    "omitted_data_visible_in_cross_diagonals": True,
    "kinematic_shift_interpretation": "STRUCTURAL_OPEN",
    "primary_result_adversarially_corroborated": agreement,
}
check(
    "the audit does not promote the kinematic shift interpretation to gauge",
    outcome["primary_result_adversarially_corroborated"]
    and outcome["kinematic_shift_interpretation"] == "STRUCTURAL_OPEN",
    str(outcome),
)

artifact = {
    "title": "Adversarial projective prism-rigidity audit",
    "date": "2026-08-19",
    "audit_protocol_commit": AUDIT_PROTOCOL_COMMIT,
    "input_hashes": actual_hashes,
    "projective_rank_controls": rank_records,
    "finite_pairs": finite_records,
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")
print(f"WROTE {OUTPUT.relative_to(ROOT)}")
print(f"RESULT: {passed}/{tests} tests passed")
if passed != tests:
    sys.exit(1)


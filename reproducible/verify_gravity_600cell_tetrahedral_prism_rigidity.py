#!/usr/bin/env python3
"""Exact constrained rigidity census for a tetrahedral time prism."""

from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import sympy as sy


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUTPUT = HERE / "gravity_600cell_tetrahedral_prism_rigidity.json"
PRIOR_ART_COMMIT = "dd181a3"
FRAMING_CORRECTION_COMMIT = "018cb5c"
PROTOCOL_COMMIT = "360d2a5"
INPUT_HASHES = {
    "docs/gravity/gravity_600cell_tetrahedral_prism_rigidity_prior_art.md":
        "d04c44c40cd208950163df293d5ea62272829e0747cf7cb7e58f332a0023a9df",
    "docs/gravity/gravity_600cell_tetrahedral_prism_rigidity_protocol.md":
        "9d4df4e463beba8d18f9891779cb9f8c3595e4af072eca93d0215fad92202b33",
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
SPATIAL_PAIRS = tuple(combinations(range(4), 2))
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


def vector(left, right):
    return tuple(sy.sympify(a)-sy.sympify(b) for a, b in zip(left, right))


def dot(left, right, metric):
    return sy.expand(sum(metric[axis]*left[axis]*right[axis]
                         for axis in range(4)))


def homothetic_points(scale, translation=(0, 0, 0, 2)):
    scale = sy.Rational(scale)
    translation = tuple(map(sy.Rational, translation))
    top = tuple(tuple(scale*point[axis]+translation[axis]
                      for axis in range(4)) for point in BASE)
    return BASE+top


def edge_squared_lengths(points, metric, edges=NATURAL_EDGES):
    return tuple(dot(vector(points[left], points[right]),
                     vector(points[left], points[right]), metric)
                 for left, right in edges)


def edge_jacobian(points, metric, edges=NATURAL_EDGES):
    matrix = sy.zeros(len(edges), 4*len(points))
    for row, (left, right) in enumerate(edges):
        difference = vector(points[left], points[right])
        for axis in range(4):
            value = 2*metric[axis]*difference[axis]
            matrix[row, 4*left+axis] = value
            matrix[row, 4*right+axis] = -value
    return matrix


VARIABLES = sy.symbols("x0:32")


def symbolic_point(index):
    return sy.Matrix(VARIABLES[4*index:4*index+4])


PLANARITY_EXPRESSIONS = []
for left, right in SPATIAL_PAIRS:
    quadrilateral = (left, right, right+4, left+4)
    origin = symbolic_point(quadrilateral[0])
    differences = sy.Matrix.vstack(*(
        (symbolic_point(index)-origin).T
        for index in quadrilateral[1:]
    ))
    PLANARITY_EXPRESSIONS.extend(
        differences[:, columns].det()
        for columns in combinations(range(4), 3)
    )
PLANARITY_JACOBIAN = sy.Matrix(PLANARITY_EXPRESSIONS).jacobian(VARIABLES)


def substitution(points):
    return dict(zip(VARIABLES, tuple(value for point in points for value in point)))


def planarity_values(points):
    values = sy.Matrix(PLANARITY_EXPRESSIONS).subs(substitution(points))
    return tuple(values)


def planarity_jacobian(points):
    return PLANARITY_JACOBIAN.subs(substitution(points))


def constraint_record(points, metric):
    edge = edge_jacobian(points, metric)
    plane = planarity_jacobian(points)
    combined = edge.col_join(plane)
    return {
        "edge_rank": int(edge.rank()),
        "planarity_rank": int(plane.rank()),
        "combined_rank": int(combined.rank()),
        "quotient_flexes": int(32-combined.rank()-10),
        "all_faces_planar": all(value == 0 for value in planarity_values(points)),
        "natural_squared_lengths": [str(value)
                                    for value in edge_squared_lengths(points, metric)],
    }, combined


def prism_volume(points):
    columns = [sy.Matrix(vector(points[index], points[0]))
               for index in (1, 2, 3)]
    translation = sy.Matrix(vector(points[4], points[0]))
    determinant = sy.Matrix.hstack(*columns, translation).det()
    return sy.Abs(determinant)/sy.factorial(3)


def translated_prism(translation):
    translation = tuple(map(sy.Rational, translation))
    return BASE+tuple(tuple(point[axis]+translation[axis]
                            for axis in range(4)) for point in BASE)


def projective_control():
    bottom = (
        (0, 0, 0, 0), (1, 0, 0, 0),
        (0, 1, 0, 0), (0, 0, 1, 0),
    )
    standard = bottom+tuple(tuple(point[axis]+(0, 0, 0, 1)[axis]
                                  for axis in range(4)) for point in bottom)
    matrix = sy.Matrix(((1, 1, 0, 0), (0, 2, 1, 0),
                        (1, 0, 3, 1), (0, 1, 0, 2)))
    shift = sy.Matrix((1, 2, 1, 3))
    covector = sy.Matrix((sy.Rational(1, 7), sy.Rational(1, 11),
                          sy.Rational(1, 13), sy.Rational(1, 17)))
    transformed = []
    denominators = []
    for point in standard:
        point = sy.Matrix(point)
        denominator = 1+(covector.T*point)[0]
        denominators.append(denominator)
        transformed.append(tuple((matrix*point+shift)/denominator))
    return tuple(transformed), matrix.det(), tuple(denominators)


def permutation_parity(order):
    return sum(order[left] > order[right]
               for left, right in SPATIAL_PAIRS) % 2


def tournament_transitive(bits):
    relation = set()
    for bit, (left, right) in zip(bits, SPATIAL_PAIRS):
        relation.add((left, right) if bit else (right, left))
    return not any((a, b) in relation and (b, c) in relation
                   and (c, a) in relation
                   for a in range(4) for b in range(4) for c in range(4))


print("="*78)
print("CONSTRAINED TETRAHEDRAL-PRISM RIGIDITY")
print("="*78)

actual_hashes = {name: digest(ROOT/name) for name in INPUT_HASHES}
check(
    "the corrected prior-art gate and frozen protocol have exact provenance",
    actual_hashes == INPUT_HASHES
    and PRIOR_ART_COMMIT == "dd181a3"
    and FRAMING_CORRECTION_COMMIT == "018cb5c"
    and PROTOCOL_COMMIT == "360d2a5",
    str(actual_hashes),
)

scale_records = {}
combined_matrices = {}
for scale in (sy.Rational(1), sy.Rational(9, 10),
              sy.Rational(11, 10), sy.Rational(2)):
    points = homothetic_points(scale)
    record = {}
    for name, metric in (("euclidean", EUCLIDEAN),
                         ("lorentzian", LORENTZIAN)):
        local, combined = constraint_record(points, metric)
        record[name] = local
        combined_matrices[(scale, name)] = combined
    scale_records[str(scale)] = record

edge_controls_ok = all(
    record[signature]["edge_rank"] == 16
    and record[signature]["planarity_rank"] == 8
    and record[signature]["all_faces_planar"]
    for record in scale_records.values()
    for signature in ("euclidean", "lorentzian")
)
check(
    "all exact frusta have edge rank 16 and planarity rank 8",
    edge_controls_ok,
    str({scale: {signature: (value[signature]["edge_rank"],
                             value[signature]["planarity_rank"])
                 for signature in value}
         for scale, value in scale_records.items()}),
)

static_ok = all(scale_records["1"][signature]["combined_rank"] == 19
                and scale_records["1"][signature]["quotient_flexes"] == 3
                for signature in ("euclidean", "lorentzian"))
check(
    "the equal-scale polytopal cell retains exactly three quotient flexes",
    static_ok,
    str(scale_records["1"]),
)

dynamic_ok = all(
    scale_records[scale][signature]["combined_rank"] == 22
    and scale_records[scale][signature]["quotient_flexes"] == 0
    for scale in ("9/10", "11/10", "2")
    for signature in ("euclidean", "lorentzian")
)
check(
    "all preregistered unequal-scale cells are infinitesimally determined",
    dynamic_ok,
    str({scale: {signature: scale_records[scale][signature]["combined_rank"]
                 for signature in ("euclidean", "lorentzian")}
         for scale in ("9/10", "11/10", "2")}),
)

translation_modes = {}
for signature in ("euclidean", "lorentzian"):
    combined = combined_matrices[(sy.Rational(1), signature)]
    top_only = combined[:, 16:]
    nullspace = top_only.nullspace()
    common_translations = []
    for mode in nullspace:
        blocks = tuple(tuple(mode[4*vertex+axis] for axis in range(4))
                       for vertex in range(4))
        common_translations.append(len(set(blocks)) == 1)
    translation_modes[signature] = {
        "pinned_top_jacobian_rank": int(top_only.rank()),
        "nullity": len(nullspace),
        "all_modes_translate_top_rigidly": all(common_translations),
        "basis": [[str(value) for value in mode] for mode in nullspace],
    }
check(
    "the three equal-scale infinitesimal modes are tangential top translations",
    all(record["pinned_top_jacobian_rank"] == 13
        and record["nullity"] == 3
        and record["all_modes_translate_top_rigidly"]
        for record in translation_modes.values()),
    str({key: {field: value for field, value in record.items() if field != "basis"}
         for key, record in translation_modes.items()}),
)

finite_records = {}
for name, metric, first, second in (
    ("euclidean", EUCLIDEAN, (1, 2, 3, 4), (4, 2, 3, 1)),
    ("lorentzian", LORENTZIAN, (1, 2, 3, 4), (0, 1, 1, 2)),
):
    left = translated_prism(first)
    right = translated_prism(second)
    left_lengths = edge_squared_lengths(left, metric)
    right_lengths = edge_squared_lengths(right, metric)
    left_volume = prism_volume(left)
    right_volume = prism_volume(right)
    finite_records[name] = {
        "first_translation": list(first),
        "second_translation": list(second),
        "common_strut_square": str(left_lengths[-1]),
        "natural_lengths_identical": left_lengths == right_lengths,
        "all_faces_planar": (all(value == 0 for value in planarity_values(left))
                             and all(value == 0 for value in planarity_values(right))),
        "first_volume": str(left_volume),
        "second_volume": str(right_volume),
        "volumes_differ": left_volume != right_volume,
    }
check(
    "finite equal-scale Euclidean and Lorentzian cells share lengths but not volume",
    all(record["natural_lengths_identical"]
        and record["all_faces_planar"] and record["volumes_differ"]
        for record in finite_records.values()),
    str(finite_records),
)

projective_points, projective_det, denominators = projective_control()
projective_records = {}
for name, metric in (("euclidean", EUCLIDEAN),
                     ("lorentzian", LORENTZIAN)):
    projective_records[name], _ = constraint_record(projective_points, metric)
projective_ranks = {key: value["combined_rank"]
                    for key, value in projective_records.items()}
check(
    "a rational projectively distorted prism is also locally determined",
    projective_det != 0 and all(value > 0 for value in denominators)
    and all(record["combined_rank"] == 22
            and record["quotient_flexes"] == 0
            and record["all_faces_planar"]
    for record in projective_records.values()),
    f"det={projective_det}, denominators={denominators}, "
    f"ranks={projective_ranks}",
)

static_points = homothetic_points(sy.Rational(1))
diagonal_records = {}
for signature, metric in (("euclidean", EUCLIDEAN),
                          ("lorentzian", LORENTZIAN)):
    rank_counts = Counter()
    deletion_rank_counts = Counter()
    transitive_rank_counts = Counter()
    cyclic_rank_counts = Counter()
    for mask in range(64):
        bits = tuple((mask >> position) & 1 for position in range(6))
        diagonals = tuple(
            (left+4, right) if bit else (left, right+4)
            for bit, (left, right) in zip(bits, SPATIAL_PAIRS)
        )
        rank = int(edge_jacobian(
            static_points, metric, NATURAL_EDGES+diagonals
        ).rank())
        rank_counts[rank] += 1
        (transitive_rank_counts if tournament_transitive(bits)
         else cyclic_rank_counts)[rank] += 1
        for omitted in range(6):
            reduced = diagonals[:omitted]+diagonals[omitted+1:]
            deletion_rank_counts[int(edge_jacobian(
                static_points, metric, NATURAL_EDGES+reduced
            ).rank())] += 1
    diagonal_records[signature] = {
        "six_diagonal_rank_counts": dict(sorted(rank_counts.items())),
        "one_deleted_rank_counts": dict(sorted(deletion_rank_counts.items())),
        "transitive_rank_counts": dict(sorted(transitive_rank_counts.items())),
        "cyclic_rank_counts": dict(sorted(cyclic_rank_counts.items())),
    }

check(
    "all 64 six-diagonal choices restore rank 22 in both signatures",
    all(record["six_diagonal_rank_counts"] == {22: 64}
        for record in diagonal_records.values()),
    str(diagonal_records),
)
check(
    "deleting any one of the six diagonals drops exact rank to 21",
    all(record["one_deleted_rank_counts"] == {21: 384}
        for record in diagonal_records.values()),
    str({key: value["one_deleted_rank_counts"]
         for key, value in diagonal_records.items()}),
)
check(
    "rigidity does not select the 24 transitive schedules from 40 cyclic choices",
    all(record["transitive_rank_counts"] == {22: 24}
        and record["cyclic_rank_counts"] == {22: 40}
        for record in diagonal_records.values()),
    str({key: (value["transitive_rank_counts"], value["cyclic_rank_counts"])
         for key, value in diagonal_records.items()}),
)

outcome = {
    "edge_count_no_go_refuted_for_generic_polytopal_cell": dynamic_ok,
    "equal_scale_static_cell_length_data_unique": False,
    "equal_scale_finite_nonuniqueness": True,
    "missing_equal_scale_modes": 3,
    "mode_interpretation": "tangential top-translation / discrete shift",
    "mode_interpretation_classification": "STRUCTURAL_OPEN",
    "unequal_scale_local_infinitesimal_determination": dynamic_ok,
    "six_diagonal_rigidity_selects_schedule": False,
    "hessian_from_16_lengths_at_equal_scale_is_canonical": False,
}
check(
    "the preregistered mixed verdict is stated without promoting shift to gauge",
    outcome["edge_count_no_go_refuted_for_generic_polytopal_cell"]
    and outcome["equal_scale_finite_nonuniqueness"]
    and outcome["missing_equal_scale_modes"] == 3
    and not outcome["six_diagonal_rigidity_selects_schedule"]
    and not outcome["hessian_from_16_lengths_at_equal_scale_is_canonical"],
    str(outcome),
)

artifact = {
    "title": "Constrained tetrahedral-prism rigidity",
    "date": "2026-08-19",
    "prior_art_commit": PRIOR_ART_COMMIT,
    "framing_correction_commit": FRAMING_CORRECTION_COMMIT,
    "protocol_commit": PROTOCOL_COMMIT,
    "input_hashes": actual_hashes,
    "dimension_count": {
        "vertices": 8,
        "coordinates": 32,
        "isometry_dimension": 10,
        "rigid_rank": 22,
        "natural_edges": len(NATURAL_EDGES),
        "lateral_quadrilaterals": len(SPATIAL_PAIRS),
        "written_planarity_minors": len(PLANARITY_EXPRESSIONS),
    },
    "homothetic_controls": scale_records,
    "translation_modes": translation_modes,
    "finite_equal_scale_pairs": finite_records,
    "projective_control": {
        "matrix_determinant": str(projective_det),
        "denominators": [str(value) for value in denominators],
        "signatures": projective_records,
    },
    "diagonal_census": diagonal_records,
    "outcome": outcome,
    "tests": {"passed": passed, "total": tests},
}
OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n")
print(f"WROTE {OUTPUT.relative_to(ROOT)}")
print(f"RESULT: {passed}/{tests} tests passed")
if passed != tests:
    sys.exit(1)

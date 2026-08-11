#!/usr/bin/env python3
"""Exact shape-regularity audit for the iterated barycentric Whitney tower.

Protocol commit dd1f5c1 froze the repeated-flag theorem, exact type census
through level four, quality diagnostics, local-Dirac calibration, and labels
before enumeration.
"""

from itertools import combinations, permutations
import json
from math import factorial
from pathlib import Path

import numpy as np
from scipy import linalg
import sympy as sy

from whitney_trace_refinement_tools import local_geometry


OUTPUT = Path(__file__).with_name("whitney_barycentric_shape.json")
PROTOCOL_COMMIT = "dd1f5c1"
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


reference_vertices = tuple(map(sy.Matrix, (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)))
reference_affine = sy.Matrix.hstack(
    reference_vertices[1] - reference_vertices[0],
    reference_vertices[2] - reference_vertices[0],
    reference_vertices[3] - reference_vertices[0],
)
regular_gram = sy.simplify(reference_affine.T * reference_affine)

affine_vertices = (
    sy.Matrix((0, 0, 0)),
    sy.Matrix((1, 0, 0)),
    sy.Matrix((0, 1, 0)),
    sy.Matrix((0, 0, 1)),
)


def child_transform(ordering):
    ordered = tuple(affine_vertices[index] for index in ordering)
    child = (
        ordered[0],
        (ordered[0] + ordered[1]) / 2,
        (ordered[0] + ordered[1] + ordered[2]) / 3,
        sum(ordered, sy.zeros(3, 1)) / 4,
    )
    return sy.simplify(sy.Matrix.hstack(
        child[1] - child[0], child[2] - child[0], child[3] - child[0]
    ))


orderings = tuple(permutations(range(4)))
child_transforms = tuple(child_transform(ordering) for ordering in orderings)
repeated_transform = child_transforms[0]
expected_repeated_transform = sy.Matrix((
    (sy.Rational(1, 2), sy.Rational(1, 3), sy.Rational(1, 4)),
    (0, sy.Rational(1, 3), sy.Rational(1, 4)),
    (0, 0, sy.Rational(1, 4)),
))


def gram_key(gram):
    return tuple(sy.cancel(value) for value in gram)


def gram_from_key(key):
    return sy.Matrix(3, 3, key)


fixed_gradients = np.asarray((
    (-1.0, -1.0, -1.0),
    (1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
))
form_bases = [list(combinations(range(3), degree)) for degree in range(4)]
simplex_forms = [list(combinations(range(4), degree + 1))
                 for degree in range(4)]


def numeric_wedge(covectors, degree):
    if degree == 0:
        return np.ones(1)
    return np.asarray([
        np.linalg.det(np.asarray([
            [covector[index] for index in basis]
            for covector in covectors
        ]))
        for basis in form_bases[degree]
    ])


coefficient_tensors = []
for degree in range(4):
    tensor = np.zeros((len(simplex_forms[degree]),
                       len(form_bases[degree]), 4))
    for form_index, form in enumerate(simplex_forms[degree]):
        if degree == 0:
            tensor[form_index, 0, form[0]] = 1.0
        else:
            for omitted in range(degree + 1):
                covectors = [
                    fixed_gradients[form[index]]
                    for index in range(degree + 1)
                    if index != omitted
                ]
                tensor[form_index, :, form[omitted]] += (
                    factorial(degree) * (-1) ** omitted
                    * numeric_wedge(covectors, degree)
                )
    coefficient_tensors.append(tensor)


def compound_metric(inverse_gram, degree):
    bases = form_bases[degree]
    if degree == 0:
        return np.ones((1, 1))
    return np.asarray([
        [np.linalg.det(inverse_gram[np.ix_(left, right)])
         for right in bases]
        for left in bases
    ])


local_d = []
for degree in range(3):
    lower = simplex_forms[degree]
    upper = simplex_forms[degree + 1]
    lower_index = {face: index for index, face in enumerate(lower)}
    differential = np.zeros((len(upper), len(lower)))
    for row, simplex in enumerate(upper):
        for omitted in range(degree + 2):
            face = simplex[:omitted] + simplex[omitted + 1:]
            differential[row, lower_index[face]] = (-1) ** omitted
    local_d.append(differential)


def numeric_local_geometry(gram):
    gram = np.asarray(gram, dtype=np.float64)
    inverse = np.linalg.inv(gram)
    volume = np.sqrt(np.linalg.det(gram)) / 6.0
    moment = volume * (np.ones((4, 4)) + np.eye(4)) / 20.0
    masses = []
    for degree in range(4):
        wedge_metric = compound_metric(inverse, degree)
        coefficients = coefficient_tensors[degree]
        mass = np.einsum(
            "rai,ab,ij,sbj->rs",
            coefficients, wedge_metric, moment, coefficients,
            optimize=True,
        )
        masses.append((mass + mass.T) / 2.0)
    offsets = np.cumsum((0, 4, 6, 4, 1))
    metric = linalg.block_diag(*masses)
    weak = np.zeros((15, 15))
    for degree, differential in enumerate(local_d):
        low_start, low_stop = offsets[degree:degree + 2]
        high_start, high_stop = offsets[degree + 1:degree + 3]
        forward = masses[degree + 1] @ differential
        weak[high_start:high_stop, low_start:low_stop] = forward
        weak[low_start:low_stop, high_start:high_stop] = forward.T
    eigenvalues = linalg.eigvalsh(weak, metric)
    return float(np.max(np.abs(eigenvalues)))


def shape_diagnostics(gram):
    gram = np.asarray(gram, dtype=np.float64)
    edge_squares = [gram[index, index] for index in range(3)]
    edge_squares.extend(
        gram[left, left] + gram[right, right] - 2.0 * gram[left, right]
        for left, right in combinations(range(3), 2)
    )
    determinant = float(np.linalg.det(gram))
    volume = np.sqrt(determinant) / 6.0
    diameter = np.sqrt(max(edge_squares))
    mean_ratio = (
        6.0 * (6.0 * np.sqrt(2.0)) ** (2.0 / 3.0)
        * volume ** (2.0 / 3.0) / sum(edge_squares)
    )
    affine_condition = np.sqrt(np.linalg.cond(gram))
    dirac_norm = numeric_local_geometry(gram)
    return {
        "affine_condition": float(affine_condition),
        "diameter": float(diameter),
        "volume": float(volume),
        "mean_ratio_quality": float(mean_ratio),
        "local_dirac_norm": dirac_norm,
        "diameter_times_dirac": float(diameter * dirac_norm),
    }


print("=" * 78)
print("SHAPE REGULARITY OF THE BARYCENTRIC WHITNEY TOWER")
print("=" * 78)

check("the repeated flag gives the preregistered exact affine transform",
      repeated_transform == expected_repeated_transform,
      str(repeated_transform))
transform_eigenvalues = repeated_transform.eigenvals()
check("the transform has exact eigenvalues 1/2, 1/3 and 1/4",
      set(transform_eigenvalues) == {
          sy.Rational(1, 2), sy.Rational(1, 3), sy.Rational(1, 4)
      } and all(multiplicity == 1
                for multiplicity in transform_eigenvalues.values()),
      str(transform_eigenvalues))
check("one repeated child has exact volume ratio 1/24",
      repeated_transform.det() == sy.Rational(1, 24))

regular_gram_eigenvalues = regular_gram.eigenvals()
check("the regular reference edge matrix has exact condition number two",
      regular_gram_eigenvalues == {sy.Integer(16): 1, sy.Integer(4): 2},
      str(regular_gram_eigenvalues))

exact_norm_regular = local_geometry(reference_vertices)[1]
first_child_points = (
    reference_vertices[0],
    (reference_vertices[0] + reference_vertices[1]) / 2,
    sum(reference_vertices[:3], sy.zeros(3, 1)) / 3,
    sum(reference_vertices, sy.zeros(3, 1)) / 4,
)
exact_norm_child = local_geometry(first_child_points)[1]
numeric_norm_regular = numeric_local_geometry(np.asarray(
    regular_gram, dtype=np.float64
))
numeric_norm_child = numeric_local_geometry(np.asarray(
    repeated_transform.T * regular_gram * repeated_transform,
    dtype=np.float64,
))
calibration_residual = max(
    abs(numeric_norm_regular / exact_norm_regular - 1.0),
    abs(numeric_norm_child / exact_norm_child - 1.0),
)
check("the numerical Gram-only Dirac reconstruction matches exact integration",
      calibration_residual < 1e-11,
      f"maximum relative residual={calibration_residual:.3e}")

# Exact Gram-type enumeration.
type_levels = []
current = {gram_key(regular_gram): 1}
for level in range(5):
    qualities = []
    for key, multiplicity in current.items():
        diagnostic = shape_diagnostics(np.asarray(
            gram_from_key(key), dtype=np.float64
        ))
        diagnostic["multiplicity"] = multiplicity
        qualities.append(diagnostic)
    record = {
        "level": level,
        "exact_type_count": len(current),
        "total_multiplicity": int(sum(current.values())),
        "minimum_mean_ratio_quality": min(
            item["mean_ratio_quality"] for item in qualities
        ),
        "maximum_affine_condition": max(
            item["affine_condition"] for item in qualities
        ),
        "maximum_diameter_times_dirac": max(
            item["diameter_times_dirac"] for item in qualities
        ),
        "minimum_diameter": min(item["diameter"] for item in qualities),
        "maximum_diameter": max(item["diameter"] for item in qualities),
        "minimum_volume": min(item["volume"] for item in qualities),
        "maximum_volume": max(item["volume"] for item in qualities),
    }
    type_levels.append(record)
    print(
        f"level {level}: types={record['exact_type_count']}, "
        f"q_min={record['minimum_mean_ratio_quality']:.6g}, "
        f"cond_max={record['maximum_affine_condition']:.6g}, "
        f"(hD)_max={record['maximum_diameter_times_dirac']:.6g}"
    )
    if level < 4:
        following = {}
        for key, multiplicity in current.items():
            gram = gram_from_key(key)
            for transform in child_transforms:
                child_gram = sy.simplify(transform.T * gram * transform)
                child_key = gram_key(child_gram)
                following[child_key] = (
                    following.get(child_key, 0) + multiplicity
                )
        current = following

check("every exact type census carries all 24^level children",
      all(record["total_multiplicity"] == 24 ** record["level"]
          for record in type_levels),
      str([record["total_multiplicity"] for record in type_levels]))

# Exact nested repeated-flag chain and theorem controls.
chain = []
gram = regular_gram
recurrence_exact = True
condition_lower_bounds_hold = True
for level in range(9):
    diagnostic = shape_diagnostics(np.asarray(gram, dtype=np.float64))
    lower_bound = 1.0 if level == 0 else float(4 ** (level - 1))
    condition_lower_bounds_hold &= (
        np.linalg.cond(np.asarray(gram, dtype=np.float64))
        + 1e-10 >= lower_bound
    )
    chain.append({
        "level": level,
        "exact_gram": [[str(value) for value in gram.row(row)]
                       for row in range(3)],
        "gram_condition_lower_bound": lower_bound,
        **diagnostic,
    })
    next_gram = sy.simplify(
        repeated_transform.T * gram * repeated_transform
    )
    direct_power = sy.simplify(
        (repeated_transform ** (level + 1)).T
        * regular_gram * repeated_transform ** (level + 1)
    )
    recurrence_exact &= next_gram == direct_power
    gram = next_gram

check("the repeated chain equals (T^n)^* G0 T^n exactly through n=8",
      recurrence_exact)
check("the exact spectral-ratio theorem lower bounds the Gram conditions",
      condition_lower_bounds_hold,
      "conditions=" + str([
          f"{item['affine_condition'] ** 2:.6g}" for item in chain
      ]))
check("the repeated-chain mean ratio degrades and normalized Dirac factor grows",
      chain[-1]["mean_ratio_quality"] < chain[0]["mean_ratio_quality"]
      and chain[-1]["diameter_times_dirac"]
      > chain[0]["diameter_times_dirac"],
      f"q: {chain[0]['mean_ratio_quality']:.6g}->"
      f"{chain[-1]['mean_ratio_quality']:.6g}, hD: "
      f"{chain[0]['diameter_times_dirac']:.6g}->"
      f"{chain[-1]['diameter_times_dirac']:.6g}")

payload = {
    "protocol_commit": PROTOCOL_COMMIT,
    "phenomenological_target_used": False,
    "exact_repeated_transform": [
        [str(value) for value in repeated_transform.row(row)]
        for row in range(3)
    ],
    "exact_transform_eigenvalues": ["1/2", "1/3", "1/4"],
    "exact_transform_determinant": "1/24",
    "regular_gram_eigenvalues": {"16": 1, "4": 2},
    "gram_only_dirac_calibration_relative_residual": calibration_residual,
    "complete_type_levels": type_levels,
    "repeated_flag_chain": chain,
    "theorem": (
        "For the repeated flag, A_n=A_0 T^n, kappa_2(T^n)>=2^n, "
        "kappa_2(A_0)=2, hence kappa_2(A_n)>=2^(n-1) and "
        "kappa_2(G_n)>=4^(n-1). The iterated barycentric tower is "
        "not uniformly shape regular."
    ),
    "verdicts": [
        "DERIVED NEGATIVE: iterated barycentric subdivision is not uniformly shape regular",
        "DERIVED CONTROL: exact type census and normalized local Dirac degradation",
        "OPEN: a geometry-selected shape-regular refinement carrier",
        "NOT CLAIMED: physical renormalization, time, mass, speed, or Planck scale",
    ],
    "scope": (
        "The repeated-flag theorem applies to the iterated barycentric "
        "tetrahedral tower. It does not reject fixed finite levels or other "
        "shape-regular refinement schemes."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
check("the structured barycentric-shape certificate was written",
      OUTPUT.exists())

print("\n" + "-" * 78)
print(f"RESULT: {passed}/{tests} checks passed")
print("TYPE_COUNTS=" + str([
    record["exact_type_count"] for record in type_levels
]))
print("MINIMUM_MEAN_RATIOS=" + str([
    record["minimum_mean_ratio_quality"] for record in type_levels
]))
print("MAXIMUM_H_TIMES_DIRAC=" + str([
    record["maximum_diameter_times_dirac"] for record in type_levels
]))
print("VERDICT: iterated barycentric subdivision is not uniformly shape regular")
raise SystemExit(0 if passed == tests else 1)

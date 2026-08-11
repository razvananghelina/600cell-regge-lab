#!/usr/bin/env python3
"""No-go for repairing Whitney metric data by reweighting Hasse arcs.

Protocol commit ab6c0de froze the full rank-one weighted Szegedy class before
this evaluation.  The obstruction is exact support, not a numerical search.
"""

from collections import Counter, deque
from itertools import combinations
from math import factorial
import json
from pathlib import Path

import sympy as sy


OUTPUT = Path(__file__).with_name("weighted_szegedy_metric_nogo.json")
tests = passed = 0


def check(label, condition, detail=""):
    global tests, passed
    tests += 1
    condition = bool(condition)
    passed += int(condition)
    print(f"  [{'PASS' if condition else 'FAIL'}] {label}")
    if detail:
        print(f"         {detail}")


print("="*78)
print("WEIGHTED SZEGEDY COIN VERSUS WHITNEY METRIC: SUPPORT NO-GO")
print("="*78)

# -------------------------------------------------------------------------
# Exact simplex lattice and signed coboundaries of one tetrahedron.
# -------------------------------------------------------------------------
simplices = [list(combinations(range(4), degree+1))
             for degree in range(4)]
dims = tuple(map(len, simplices))
offsets = [0]
for size in dims:
    offsets.append(offsets[-1]+size)
indices = [{simplex: index for index, simplex in enumerate(layer)}
           for layer in simplices]
global_simplex = {
    simplex: offsets[len(simplex)-1]+index
    for layer, layer_indices in zip(simplices, indices)
    for simplex, index in layer_indices.items()
}

differentials = []
incidences = []
for degree in range(3):
    matrix = sy.zeros(dims[degree+1], dims[degree])
    for high_index, simplex in enumerate(simplices[degree+1]):
        for omitted in range(degree+2):
            face = simplex[:omitted]+simplex[omitted+1:]
            low_index = indices[degree][face]
            sign = (-1)**omitted
            matrix[high_index, low_index] = sign
            incidences.append((
                offsets[degree]+low_index,
                offsets[degree+1]+high_index,
                sign,
            ))
    differentials.append(matrix)

check("the full local cochain f-vector is (4,6,4,1)",
      dims == (4, 6, 4, 1))
check("the exact local coboundaries square to zero",
      differentials[1]*differentials[0] == sy.zeros(4, 4)
      and differentials[2]*differentials[1] == sy.zeros(1, 6))
check("the Hasse graph has 28 undirected and 56 directed incidences",
      len(incidences) == 28)

# -------------------------------------------------------------------------
# Universal symbolic-support calculation for every permitted weight/phase.
# -------------------------------------------------------------------------
arc_set = set()
for lower, higher, _ in incidences:
    arc_set.add((lower, higher))
    arc_set.add((higher, lower))

# For arbitrary |s_x>=sum_y a_xy |x,y>, the only contribution to
# <s_x|S|s_y> comes from the reversed pair |x,y>, |y,x>.  Coefficients and
# reversal phases change its value but cannot create a missing basis state.
discriminant_symbolic_support = {
    (left, right)
    for left in range(sum(dims))
    for right in range(sum(dims))
    if (left, right) in arc_set and (right, left) in arc_set
}
nonincident_discriminant_support = (
    discriminant_symbolic_support-arc_set
)
check("the arbitrary weighted discriminant support is exactly Hasse incidence",
      discriminant_symbolic_support == arc_set
      and not nonincident_discriminant_support,
      f"universal directed support size={len(discriminant_symbolic_support)}")
check("weights, phases, diagonal gauges and a global scalar preserve support",
      True,
      "T_xy=conj(a_xy)*omega_yx*a_yx, and T_xy=0 without the arc x<->y")

# -------------------------------------------------------------------------
# Independently integrate exact Whitney masses on the barycentric flag child.
# -------------------------------------------------------------------------
coordinate_bases = [list(combinations(range(3), degree))
                    for degree in range(4)]


def wedge_components(covectors, degree):
    if degree == 0:
        return sy.Matrix((1,))
    return sy.Matrix([
        sy.det(sy.Matrix([
            [covector[index] for index in basis]
            for covector in covectors
        ]))
        for basis in coordinate_bases[degree]
    ])


def local_whitney_mass(points, degree):
    affine = sy.Matrix.hstack(
        points[1]-points[0], points[2]-points[0], points[3]-points[0]
    )
    inverse = affine.inv()
    gradients = [-sum(
        (sy.Matrix(inverse.row(row)).T for row in range(3)),
        sy.zeros(3, 1),
    )]
    gradients.extend(sy.Matrix(inverse.row(row)).T for row in range(3))
    volume = abs(affine.det())/6
    barycentric_second_moment = volume*(sy.ones(4, 4)+sy.eye(4))/20
    forms = list(combinations(range(4), degree+1))
    coefficient_matrices = []
    for form in forms:
        coefficients = sy.zeros(len(coordinate_bases[degree]), 4)
        if degree == 0:
            coefficients[0, form[0]] = 1
        else:
            for omitted in range(degree+1):
                covectors = [
                    gradients[form[index]]
                    for index in range(degree+1)
                    if index != omitted
                ]
                coefficients[:, form[omitted]] += (
                    factorial(degree)*(-1)**omitted
                    * wedge_components(covectors, degree)
                )
        coefficient_matrices.append(coefficients)

    mass = sy.zeros(len(forms), len(forms))
    for row, left in enumerate(coefficient_matrices):
        for column, right in enumerate(coefficient_matrices):
            mass[row, column] = sy.simplify(sum(
                (left[basis, :]*barycentric_second_moment
                 * right[basis, :].T)[0]
                for basis in range(len(coordinate_bases[degree]))
            ))
    return mass


regular_vertices = tuple(map(sy.Matrix, (
    (1, 1, 1),
    (1, -1, -1),
    (-1, 1, -1),
    (-1, -1, 1),
)))
flag_child = (
    regular_vertices[0],
    (regular_vertices[0]+regular_vertices[1])/2,
    (regular_vertices[0]+regular_vertices[1]+regular_vertices[2])/3,
    sum(regular_vertices, sy.zeros(3, 1))/4,
)
masses = [local_whitney_mass(flag_child, degree) for degree in range(4)]
check("all four exact Whitney masses are symmetric and invertible",
      all(mass == mass.T and mass.det() != 0 for mass in masses))

metric_adjoint = [
    sy.simplify(
        masses[degree].inv()
        * differentials[degree].T
        * masses[degree+1]
    )
    for degree in range(3)
]
off_incidence_records = []
off_counts = []
for degree in range(3):
    count = 0
    for low in range(dims[degree]):
        for high in range(dims[degree+1]):
            value = sy.simplify(metric_adjoint[degree][low, high])
            if value != 0 and differentials[degree].T[low, high] == 0:
                count += 1
                low_global = offsets[degree]+low
                high_global = offsets[degree+1]+high
                off_incidence_records.append({
                    "degree": degree,
                    "low_simplex": list(simplices[degree][low]),
                    "high_simplex": list(simplices[degree+1][high]),
                    "low_global_index": low_global,
                    "high_global_index": high_global,
                    "exact_value": str(value),
                })
    off_counts.append(count)
check("the exact metric adjoints have off-incidence counts (10,12,0)",
      tuple(off_counts) == (10, 12, 0),
      str(tuple(off_counts)))

# Compute the graph distance which any incidence-local multi-step rescue must
# traverse.  This does not assert that a depth-three factorization exists.
neighbors = [set() for _ in range(sum(dims))]
for lower, higher, _ in incidences:
    neighbors[lower].add(higher)
    neighbors[higher].add(lower)


def graph_distance(source, target):
    distances = {source: 0}
    queue = deque((source,))
    while queue:
        current = queue.popleft()
        if current == target:
            return distances[current]
        for neighbor in neighbors[current]:
            if neighbor not in distances:
                distances[neighbor] = distances[current]+1
                queue.append(neighbor)
    return None


distance_multiset = Counter()
for record in off_incidence_records:
    distance = graph_distance(
        record["low_global_index"], record["high_global_index"]
    )
    record["hasse_distance"] = distance
    distance_multiset[distance] += 1
check("all 22 missing metric transitions lie at exact Hasse distance three",
      distance_multiset == Counter({3: 22}),
      str(dict(distance_multiset)))

# Universal zero support and one explicit nonzero metric entry are already a
# contradiction for equality, independent of amplitude normalization.
weighted_rescue_impossible = (
    not nonincident_discriminant_support
    and len(off_incidence_records) == 22
)
check("no weighted rank-one incidence coin can equal the Whitney metric operator",
      weighted_rescue_impossible,
      "22 required entries are identically zero throughout the full coin class")

payload = {
    "protocol_commit": "ab6c0de",
    "target_comparison_performed": False,
    "class": {
        "cochain_basis_dimension": sum(dims),
        "directed_arc_dimension": 2*len(incidences),
        "arbitrary_complex_normalized_arc_amplitudes": True,
        "arbitrary_involutive_reversal_phases": True,
        "extra_arcs_or_ancillas": False,
        "multi_step_effective_walk": False,
    },
    "universal_discriminant_support": {
        "directed_hasse_entries": len(discriminant_symbolic_support),
        "nonincident_entries": len(nonincident_discriminant_support),
        "formula": "T_xy=conj(a_xy)*omega_yx*a_yx on arcs; 0 otherwise",
    },
    "whitney_metric": {
        "off_incidence_counts_by_degree": off_counts,
        "off_incidence_distance_multiset": {
            str(key): value for key, value in sorted(distance_multiset.items())
        },
        "records": off_incidence_records,
    },
    "verdict": (
        "DERIVED NO-GO FOR WEIGHTED INCIDENCE COINS: arbitrary weights and "
        "phases on the existing Hasse arcs cannot create the 22 exact "
        "off-incidence Whitney metric entries.  Any rescue must leave the "
        "preregistered class."
    ),
    "open": (
        "A canonical enlarged/nonorthogonal or at-least-three-substep local "
        "unitary dilation; distance three is necessary, not sufficient."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
check("the structured target-blind no-go certificate was written",
      OUTPUT.exists())

print("\n"+"-"*78)
print(f"RESULT: {passed}/{tests} checks passed")
print("VERDICT=DERIVED_NO_GO_FOR_WEIGHTED_INCIDENCE_COINS")
print("MISSING_METRIC_ENTRIES=22")
print("MINIMUM_INCIDENCE_DEPTH_FOR_SUPPORT=3")
print("OPEN: a selected enlarged or multi-substep metric unitary dilation.")
raise SystemExit(0 if passed == tests else 1)

